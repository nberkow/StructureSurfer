# StructureSurfer
A tool for browsing RNA secondary structure information

Structure Surfer is a web tool for scientists who want to browse RNA secondary structure information from different labs. It's online here:

http://tesla.pcbi.upenn.edu/structuresurfer

This repository contains the database as a mySQL dump file for users who want to install it themselves rather than use the web tool. It also contains the Python script that the webtool uses to browse the database.

#### Creating The Database
You will need mySQL installed and a user account with the ability to GRANT SELECT priviledges. Run this command:
>mysql -p < structure_surfer.mysql

This makes an empty database with three tables:
>structure_score - RNA secondary structure scores with genomic coordinates

>structure_source - The experiments that generates the scores

>transcript - Exon coordinates for transcripts

#### Populating the database
# TODO

#### Browsing The Database
structurePlotMaker.py is a tool for browsing the database. It can take a few types of requests. It generates a table of results in plaintext and in xml and an xml plot.

###### Get scores from all datasets using genomic coordinate
>python2.7 makeStructurePlot.py -c chr7 -s 45459777 -e 45459811 -g mm -pfx my_output_file

-s and -e The start and end coordinates

-c The chromosome

-g Specifies the genome: mm (mouse), hs (human) or at (thale cress)

-pfx The prefix for the three output files

###### Get average score of several regions using a bed file
In some cases it's useful to take several regions of interest and find the average score profile across them. 

>python2.7 makeStructurePlot.py -b my_input_file.bed -g mm -pfx my_output_file

-b The file name of a bed file containing the coordinates of interest. All bed intervals must be of the same size. 

###### Get scores from an annotated transcript
>python2.7 makeStructurePlot.py -t AT3G61897.1 -g at -pfx my_output_file

-t Transcript ID. This must be a transcript ID that exists in the transcript table


