###################################################
## Author: Zac Brown <zac@zacbrown.org>
## Date: 01.28.2009
## File: biodb_populator.py
## Purpose: populate local BioSQL db with sequences
###################################################

from processing import Pool
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
import yaml
import os
import sys

Entrez.email = "zbrown@miami.edu"

gen_path = "genomes-all"
config_path = "config.yml"
genome_file_list = os.listdir(gen_path)

# load the config
config_file_h = open(config_path, "rU")
conf_vals = yaml.load(config_file_h)
db_vals = conf_vals["default"]

#open a db connection
server = BioSeqDatabase.open_database(driver="MySQLdb",
                                      user=db_vals["username"],
                                      passwd=db_vals["password"],
                                      host=db_vals["server"], db="biosql")

db_h = server.new_database("dp_prob",
                           description="Taxonomy relations for profiled genomes")
server.adaptor.commit()

for file in genome_file_list:
    path_h = open(gen_path+"/"+file, "rU")

    for record in SeqIO.parse(path_h, "fasta"):
        entrez_h = Entrez.efetch(id=record.id, db="nucleotide",
                                 rettype="genbank")
        print "Fetched: " + record.description
        db_h.load(SeqIO.parse(entrez_h, "genbank"))
        print "Loaded: " + record.description + "\n"

server.adaptor.commit()
