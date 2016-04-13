###############################################################
#
# StructurePlotMaker
#
# This object contains methods for querying the structure
# database. Currently supports
# - genomic coordinate queries
# - bedfile upload queries
# - gene id queries


import MySQLdb as mdb
import pygal
import sys

class StructurePlotMaker:

	
	def __init__(self):
		self.bed_intervals = []
		self.genomic_coords = []
		self.interval_len = 0

		"""
		# Deprecated code for displaying the nucleotide
		self.single_template = ("SELECT ss.pos-%s as relpos, nuc.nt, score, -1 FROM structure_score ss " + 
					"JOIN structure_source src on src.source_id = ss.source_id " + 
                                        "LEFT OUTER JOIN nucleotide nuc ON " +
                                        "src.spec_id = nuc.species_id " +
                                        "AND ss.pos = nuc.pos AND ss.chrom = nuc.chrom " +
					"WHERE src.source_id = %s and ss.chrom = '%s' " + 
					"AND ss.pos BETWEEN %s AND %s - 1")
		"""

		self.single_template = ("SELECT ss.pos-%s as relpos, 'N', score, -1 FROM structure_score ss " + 
					"JOIN structure_source src on src.source_id = ss.source_id " + 
					"WHERE src.source_id = %s and ss.chrom = '%s' " + 
					"AND ss.pos BETWEEN %s AND %s - 1")

		self.template =  ("(SELECT pos-%s as relpos, 'N', score " +
					"FROM structure_score " +
					"WHERE source_id = %s " +
					"AND chrom = '%s' AND pos " +
					"BETWEEN %s AND %s - 1)")
		self.query = ""
		self.con = mdb.connect('127.0.0.1', 'ss_user', 'ea8dd48b0a35c9c92281b68193eec33b', 'structure_surfer', unix_socket='/var/lib/mysql/mysql.sock')

		self.MAX_BED_ROWS = 100000

	## Read in a bed file and make sure all intervals are the same length
	## Create a SQL query to get the scores
	def bed_to_query(self, bed_file, source_id):

		i = 0
		query_blocks = []		
		for bed in bed_file:
			i += 1
			if i >= self.MAX_BED_ROWS:
				e = "too many rows"
				return (False, e)


			elts = bed.split("\t")
			if len(elts) < 3:
				e = "line with fewer than 3 elements:   " + bed
				return (False, e) 

			(chrm, start, end) = elts[0], int(elts[1]), int(elts[2])
			length = end - start
			if self.interval_len == 0:
				self.interval_len = length
			if self.interval_len != length:
				e = "print inconsistent length: %s v %s:    " % (length, self.interval_len) + bed
				return (False, e)

			query_blocks.append(self.template % (start, source_id, chrm, start, end))


		self.query = "\nSELECT relpos, 'N', avg(score), std(score) from (" + "\nUNION\n".join(query_blocks) + "\n) t_union group by relpos;"
		#print self.query
		return (True, None)

        def coord_to_query(self, chrm, start, end, source_id):
                self.query = self.single_template % (start, source_id, chrm, start, end)
                self.interval_len = end - start

	def transid_to_query(self, trans_id, source_id):

		self.query = """
				SELECT ss.pos - x.start as relpos, 'N', ss.score, -1 FROM structure_score ss 
				JOIN transcript x 
				ON ss.chrom = x.chrom AND ss.pos BETWEEN x.start AND x.end 
				WHERE trans_id = '%s' and ss.source_id = %s
				ORDER BY start;""" % (trans_id, source_id) 	

	def get_abs_exon_coords(self, trans_id):
		cursor = self.con.cursor()
		x_len_query = "SELECT start, end, strand FROM transcript WHERE trans_id = '%s' ORDER BY start" % (trans_id)
		cursor.execute(x_len_query)
		res = cursor.fetchall()
		
		for r in res:
			self.genomic_coords = self.genomic_coords + range(r[0], r[1])			
		self.interval_len = len(self.genomic_coords)

	def run_query(self):
		cursor = self.con.cursor()
		cursor.execute(self.query)
		res = cursor.fetchall()
		#print res
		data_set = {}
		for r in res:
                        n = 'N'
			if r[1]:
				n = r[1]	  
		
			#data_set[int(r[0])] = (n, round(r[2],3), round(r[3],3))
			data_set[int(r[0])] = (round(r[2],3), round(r[3],3))

		return data_set

	def iterate_through_matches(self, motif_id, source_id, flank, genome):
		top_cursor = self.con.cursor()
		cursor = self.con.cursor()

		alpha = 10**-10
		query = "SELECT chr, start - %s, stop + %s FROM motif_matches WHERE motif_ID = '%s' AND species = '%s' AND p_val <= %s;" 
		top_cursor.execute(query % (flank, flank, motif_id, genome, alpha))

		sums = {}
		counts = {}
		square_sums = {}
		data_set = {}
		means = {}
		result = top_cursor.fetchone()

		while result:
			chrm, start, end = result
			self.coord_to_query(chrm, int(start), int(end), source_id)
			cursor.execute(self.query)

			res = cursor.fetchone()
			while res:
				relpos, score, f = res
				if not counts.has_key(relpos):
					counts[relpos] = 0
					sums[relpos] = 0.
					square_sums[relpos] = 0.
				sums[relpos] += score
				counts[relpos] += 1
				res = cursor.fetchone()
			result = top_cursor.fetchone()

		for s in counts.keys():
			means[s] = (sums[s]/counts[s])
	
		top_cursor.execute(query % (flank, flank, motif_id, genome, alpha))
		result = top_cursor.fetchone()

		while result:
			chrm, start, end = result 
			self.coord_to_query(chrm, int(start), int(end), source_id)
			cursor.execute(self.query)
			res = cursor.fetchone()
			while res:
				relpos, score, f = res
				square_sums[relpos] += (score - means[relpos])**2
				res = cursor.fetchone()
			result = top_cursor.fetchone()	

		for s in counts.keys():
			sd = (square_sums[s]/counts[s])**0.5
			data_set[s] = (means[s], sd)		

		return data_set				

	def iterate_through_bed(self, bed_file, source_id):
		cursor = self.con.cursor()
		b = bed_file.readline()
		elts = b.split("\t")

		chrm, start, end = elts[0:3]
		width = int(end) - int(start)
		sums = {}
		counts = {}
		square_sums = {}
		data_set = {}
		means = {}

		while b:
			elts = b.split("\t")
			chrm, start, end = elts[0:3]
			w = int(end) - int(start)
			if w != width:
				print "BED error. all intervals must be of equal length: %s != %s, %s\t%s\t%s" % (w, width, chrm, start, end)
				sys.exit()
			self.coord_to_query(chrm, int(start), int(end), source_id)
			cursor.execute(self.query)
			res = cursor.fetchone()
			while res:
				relpos, n, score, f = res
				if not counts.has_key(relpos):
					counts[relpos] = 0
					sums[relpos] = 0.
					square_sums[relpos] = 0.
				sums[relpos] += score
				counts[relpos] += 1
				res = cursor.fetchone()
			b = bed_file.readline()

		for s in counts.keys():
			means[s] = (sums[s]/counts[s])
		
		bed_file.seek(0)
		b = bed_file.readline()
		while b:
			elts = b.split("\t")
			chrm, start, end = elts[0:3]
			self.coord_to_query(chrm, int(start), int(end), source_id)
			cursor.execute(self.query)
			res = cursor.fetchone()
			while res:
				relpos, n, score, f = res
				square_sums[relpos] += (score - means[relpos])**2
				res = cursor.fetchone()
			b = bed_file.readline()

		for s in counts.keys():
			sd = (square_sums[s]/counts[s])**0.5
			data_set[s] = (means[s], sd)		

		return data_set				



















