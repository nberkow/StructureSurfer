########################################
#
# makeStructurePlot
#
# this script is called by the StructureSurfer
# web tool whenever a query is submitted.
# It takes input from the web interface and creates
# a StructurePlotMaker object to query the database.
# It formats the results of the query and makes the
# output that the tool displays
#
#

import StructurePlotMaker
import MotifStructureTool
import os.path
import pygal
from pygal.style import DefaultStyle
import argparse
import sys
import re

class StructurePlotHandler:

    def __init__(self):

	# Currently available datasets
        self.structure_sources = {
           1 : "PARS_rep1",
	   2 : "PARS_rep2" ,    
	   3 : "icSHAPE_vitro",  
	   4 : "icSHAPE_vivo", 
	   5 : "DMS",   
           6 : "ds_ss_RNAse_Seq",
           7 : "atDMS_rep1",
           8 : "atDMS_rep2",
           9 : "ds_ss_RNAse_Seq",
          10 : "At_li_et_al"
          #11 : "At_polyA_2",
          #12 : "At_polyA_1"
         }

        self.invert = set([3,4,5])
        self.sources = []
        self.species = 0

    def handle_coord_request(self, chrom, start, end, genome, prefix):

        if genome == 'hs':
             self.sources = [1, 2, 5, 6]
             self.species_id = 1
        elif genome== 'mm':
             self.species_id = 2
             self.sources = [3, 4]
        elif genome == 'at':
             self.species_id = 3
             self.sources = [7, 8, 9, 10] #, 11, 12]
        else:
             print "unrecognized genome: '%s'" % (genome)
             sys.exit()

        spm = StructurePlotMaker.StructurePlotMaker()
        summ_data = self.run_coord_search(chrom, start, end, spm)
        interval_len = end - start
        summary = self.format_report(summ_data, interval_len)
        graph = self.make_summary_graph(summ_data, interval_len)
        report = self.format_for_load(summ_data, interval_len, range(start,end))

        g = open(prefix + "_graph.xml", 'w+')
        print >> g, graph
        g.close()

        s = open(prefix + "_summary.xml", 'w+')
        print >> s, summary
        s.close()

        t = open(prefix + ".tbl", 'w')
        print >> t, report
        t.close()


    def handle_bed_request(self, bed_fname, genome, prefix):

        if genome == 'hs':
             self.sources = [1, 2, 5, 6]
        elif genome== 'mm':
             self.sources = [3, 4]
        elif genome == 'at':
             self.species_id = 3
             self.sources = [7, 8, 9, 10] #, 11, 12]
        else:
             print "unrecognized genome: '%s'" % (genome)
             sys.exit()

        spm = StructurePlotMaker.StructurePlotMaker()
        summ_data = self.run_bed_search_unbatched(open(bed_fname, 'r'), spm)
        summary = self.format_report(summ_data, spm.interval_len)
        graph = self.make_summary_graph(summ_data, spm.interval_len)
	report = self.format_for_load(summ_data, spm.interval_len)

        g = open(prefix + "_graph.xml", 'w+')
        print >> g, graph
        g.close()

        s = open(prefix + "_summary.xml", 'w+')
        print >> s, summary 
        s.close()

	t = open(prefix + ".tbl", 'w')
	print >> t, report
	t.close()


    def handle_transcript_id_request(self, transcript_id, genome, prefix):

        if genome == 'hs':
             self.sources = [1, 2, 5, 6]
        elif genome== 'mm':
             self.sources = [3, 4]
        elif genome == 'at':
             self.species_id = 3
             self.sources = [7, 8, 9, 10] #, 11, 12]
        else:
             print "unrecognized genome: '%s'" % (genome)
             sys.exit()

	spm = StructurePlotMaker.StructurePlotMaker()
	summ_data = self.run_transcript_search(transcript_id, spm)
        summary = self.format_report(summ_data, spm.interval_len)
        graph = self.make_summary_graph(summ_data, spm.interval_len)
	report = self.format_for_load(summ_data, spm.interval_len, spm.genomic_coords)

        g = open(prefix + "_graph.xml", 'w+')
        print >> g, graph
        g.close()

        s = open(prefix + "_summary.xml", 'w+')
        print >> s, summary 
        s.close()

	t = open(prefix + ".tbl", 'w')
	print >> t, report
	t.close()


    # Currently depricated
    def handle_motif_request(self, motif_id, genome, prefix):

        if genome == 'hs':
             self.sources = [1, 2, 5, 6]
        elif genome== 'mm':
             self.sources = [3, 4]
        else:
             print "unrecognized genome"
             sys.exit()

        flank = 20 
        spm = StructurePlotMaker.StructurePlotMaker()
        summ_data = self.run_motif_search_unbatched(motif_id, spm, flank, genome)
        summary = self.format_report(summ_data, spm.interval_len)
        graph = self.make_summary_graph(summ_data, spm.interval_len)
	report = self.format_for_load(summ_data, spm.interval_len)

        g = open(prefix + "_graph.xml", 'w+')
        print >> g, graph
        g.close()

        s = open(prefix + "_summary.xml", 'w+')
        print >> s, summary
        s.close()

        d = open(prefix + ".tbl", 'w')
        print >> d, report 
        d.close()
 
    def run_bed_search(self, bed, spm):
        data_sets = []	

        for s in self.sources:
            bed.seek(0)
            (validBed, e) = spm.bed_to_query(bed, s)
            if not validBed:
                print "BED file error"
                print e
                sys.exit()
            else:
                res = spm.run_query()
  	        data_sets.append(res)
        bed.close()
        return data_sets

    def run_bed_search_unbatched(self, bed, spm):
        data_sets = []
	
        for s in self.sources:
            bed.seek(0)
            data_set = spm.iterate_through_bed(bed, s)
            data_sets.append(data_set)
        return data_sets

    def run_transcript_search(self, transcript_id, spm):
        data_sets = []
        for s in self.sources:
            spm.transid_to_query(transcript_id, s)
            spm.get_abs_exon_coords(transcript_id)
            res = spm.run_query()
            data_sets.append(res)
        return data_sets

    def run_motif_search_unbatched(self, motif_id, spm, flank, g):

        data_sets = []
	genomes = {'hs':'Homo_sapiens', 'mm':'Mus_musculus'}
	if not g in genomes.keys():
		print "unrecognized genome"
		sys.exit()
	
        for s in self.sources:
            data_set = spm.iterate_through_matches(motif_id, s, flank, genomes[g])
            data_sets.append(data_set)
	return data_sets

    def run_coord_search(self, chrom, start, end, spm):
        data_sets = []
        for s in self.sources:
            spm.coord_to_query(chrom, start, end, s)
            res = spm.run_query()
            data_sets.append(res)
        return data_sets


    def format_report(self, data_sets, interval_length):

       summary = '<table border = "1"><b>\n'
       summary = summary + '\t<tr><td>pos</td>\n'

       for s in self.sources:
           summary = summary + '\t\t<td>%s Avg</td><td>%s StDev</td>' % (self.structure_sources[s], self.structure_sources[s])
       summary = summary + '\t</tr></b>\n'
       sets = len(data_sets)

       for x in range(interval_length):
           summary = summary + '\t<tr><td>%s</td>\n' % (x)
           for d in range(sets):
               scores = ('NA', 'NA', 'NA')
               s0 = 'NA'
               s1 = 'NA'
               if data_sets[d].has_key(x):
                   scores = data_sets[d][x]
               if scores[0] != 'NA':
                   s0 = "%.2f" % scores[0]
               if scores[1] != 'NA':
                   s1 = "%.2f" % scores[1]
               summary = summary + '\t\t<td>%s</td><td>%s</td>' % (s0, s1) 
           summary = summary + '\t</tr>'

       summary += '</table>'
       return summary


    def format_for_load(self, data_sets, interval_length, genome_pos_list=None):

       summary = ''
       summary = summary + 'relative pos\tchromosomal pos'

       for s in self.sources:
           summary = summary + '\t%s Avg\t%s StDev' % (self.structure_sources[s], self.structure_sources[s])
       summary = summary + '\n'
       sets = len(data_sets)

       for x in range(interval_length):
           if genome_pos_list:
               g = genome_pos_list[x]
           else:
               g = 'NA'

           summary = summary + '%s\t%s' % (x, g)
           for d in range(sets):
               scores = ('NA', 'NA', 'NA')
               if data_sets[d].has_key(x):
                   scores = data_sets[d][x]
               summary = summary + '\t%s\t%s' % (scores[0], scores[1]) 
           summary = summary + '\n'

       return summary

    def organize_scores_for_plot(self, data_sets, interval_length):
        plot_averages = {}
        plot_sds = {}

        for x in range(len(data_sets)):
            source_id = self.sources[x]
            source_name = self.structure_sources[source_id]
            data_set = data_sets[x]

            plot_averages[source_name] = [None] * interval_length
            plot_sds[source_name] = [None] * interval_length
            for y in data_set.keys():
                plot_averages[source_name][y] = data_set[y][0]
                plot_sds[source_name][y] = data_set[y][1]

        return(plot_averages, plot_sds)

    def organize_and_scale_scores_for_plot(self, data_sets, interval_length):
        plot_averages = {}
        plot_sds = {}

        for x in range(len(data_sets)):
            source_id = self.sources[x]
            source_name = self.structure_sources[source_id]
            data_set = data_sets[x]

            plot_averages[source_name] = [None] * interval_length
            score_sum = 0.
            valid_scores = len(data_set.keys())

            # Organize the data into lists with missing values filled in
            plot_sds[source_name] = [None] * interval_length
            for y in data_set.keys():
                plot_averages[source_name][y] = data_set[y][0]
                plot_sds[source_name][y] = data_set[y][1]
                score_sum += data_set[y][0]

            # Calculate the scaling factor
            offset = 01
            if valid_scores > 0:
                offset = score_sum/valid_scores
            square_diff_sum = 0.
            if valid_scores > 0:
                for y in data_set.keys():
                    square_diff_sum += (offset - data_set[y][0])**2
            if square_diff_sum > 0 and valid_scores > 0:
                    scale = (square_diff_sum/valid_scores)**0.5
            else:
                    scale = 1

            print scale
            for s in range(len(plot_averages[source_name])):
                if plot_averages[source_name][s] != None and valid_scores > 0:
                    plot_averages[source_name][s] = (plot_averages[source_name][s] - offset)/scale

                    if source_id in self.invert:
                        plot_averages[source_name][s] = 0 - plot_averages[source_name][s]
        return(plot_averages, plot_sds)

    def make_summary_graph(self, data_sets, interval_length):

        #(plot_averages, plot_sds) = self.organize_scores_for_plot(data_sets, interval_length)
        (plot_averages, plot_sds) = self.organize_and_scale_scores_for_plot(data_sets, interval_length)

	line_chart = pygal.Line(style=DefaultStyle, 
		#interpolate='hermite',
		y_title = 'Relative Structure', 
		x_title = 'Position', 
		show_minor_x_labels=False, 
		show_dots=True,
		stroke=True,
		x_label_rotation=-90)
        line_chart.title = 'Average Structure Profile'
        line_chart.x_labels = map(str, range(interval_length))
        line_chart.x_labels_major = map(str, range(0, interval_length, 10))

        for x in plot_averages.keys():
             line_chart.add(x, plot_averages[x])
        #graph = line_chart.render(disable_xml_declaration=True)
        #return "<figure>" + graph + "</figure>"
	graph = line_chart.render()
	return graph


if __name__ == "__main__":

     parser = argparse.ArgumentParser(description='Pull up genomic coordinate from a bed file or the command line and return structure data')
     parser.add_argument('-c', '--chromosome')
     parser.add_argument('-g', '--genome', help="genome matching bed file or coordinates (currently supports 'hs' and 'mm'")
     parser.add_argument('-s', '--start-position', help="start position of the window in absolute chromosomal coordinates")
     parser.add_argument('-e', '--end-position', help="end position of the window in absolute chromosomal coordinates")
     parser.add_argument('-b', '--bed-file', help="a bed file of positions to determine average structure. All bed intervals must be of equal length")
     parser.add_argument('-x', '--transcript-id', help="A transcript ID for getting the structure scores of a spliced transcript")
     #parser.add_argument('-m', '--motif-id', help="motif id of an RBP. Scores will be aggregated across all motif sites ")
     parser.add_argument('-pfx', '--outfile-prefix', help="a prefix for the ouput files")
     args = parser.parse_args()

     coord_entered = (args.chromosome and args.start_position and args.end_position)
     sph = StructurePlotHandler()

     if not args.bed_file and not coord_entered and not args.transcript_id: 
         print "please enter genomic coordinates, and RBP motif ID or a bed file"
         sys.exit()

     if  args.bed_file and coord_entered :
         print "please enter coordinates OR a bed file"
         sys.exit()

     if not args.genome:
         print "please specify a genome (hs or mm)"
         sys.exit()

     if not args.outfile_prefix:
         print "please choose a prefix for the output files"
         sys.exit()

     if coord_entered:
         numeric = re.compile(r'[^\d]+')

         sph.handle_coord_request(args.chromosome, int(numeric.sub("", args.start_position)), int(numeric.sub("", args.end_position)), args.genome, args.outfile_prefix)

     if args.bed_file:
         sph.handle_bed_request(args.bed_file, args.genome, args.outfile_prefix)

     #if args.motif_id:
     #    sph.handle_motif_request(args.motif_id, args.genome, args.outfile_prefix)

     if args.transcript_id:
	sph.handle_transcript_id_request(args.transcript_id, args.genome, args.outfile_prefix)





