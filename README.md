# StructureSurfer
A tool for browsing RNA secondary structure information

Structure Surfer is a web tool for scientists who want to browse RNA secondary structure information from different labs. It's online here:

http://tesla.pcbi.upenn.edu/structuresurfer

This repository contains the database as a mySQL dump file for users who want to install it themselves rather than use the web tool. It also contains the Python script that the webtool uses to browse the database.

## Creating The Database
You will need mySQL installed and a user account with the ability to GRANT SELECT priviledges. Run this command:

mysql -p < structure_surfer.mysql

This creates a database with three tables:

structure_score 

structure_source

transcript      

