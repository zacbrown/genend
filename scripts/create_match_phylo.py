#!/usr/bin/env python

import MySQLdb
import yaml
import sys
import datetime
import re

date = datetime.date.today().isoformat()

if len(sys.argv) < 3:
    print "ERROR while doin' work son!"
    sys.exit(-1)

# load the config
config_path = sys.argv[2]
config_file_h = open(config_path, "rU")
conf_vals = yaml.load(config_file_h)
db_vals = conf_vals["default"]

def fetchTaxonomy (in_seq_name, cursor):
    cur_taxon_id = ""; next_taxon_id = ""
    ret_dict = {}

    seq_name = in_seq_name.split('.')[0].strip("\n")

    # get current taxon id
    query = "SELECT taxon_id FROM bioentry WHERE name='%s'" % seq_name

    cursor.execute(query)
    cur_taxon_id = cursor.fetchone()

    if cur_taxon_id != None: cur_taxon_id = cur_taxon_id[0]
    else: return None

    while True:
        query = "SELECT node_rank, parent_taxon_id FROM taxon WHERE taxon_id='%s'" % cur_taxon_id
        cursor.execute(query)
        (rank, next_taxon_id) = cursor.fetchone()

        if next_taxon_id == 1: break

        ret_dict[rank] = str(cur_taxon_id)

        cur_taxon_id = next_taxon_id

    return ret_dict




# load file
file_h = open(sys.argv[1], "rU")

db_h = MySQLdb.connect(host=db_vals["server"], user=db_vals["username"],
                       passwd=db_vals["password"], db=db_vals["database"])

match_vals = {"superkingdom":0, "phylum":0, "class":0, "order":0, "family":0, "genus":0}
total_vals = {"superkingdom":0, "phylum":0, "class":0, "order":0, "family":0, "genus":0}

cursor = db_h.cursor()
count_measured = 0
counter = 0
cur_kmer = 3

splitter = re.compile(r'([-.])')
piece_size = splitter.split(sys.argv[1])[2]
genus_h = open("genus-"+piece_size+"-"+date+".dat",'w')
family_h = open("family-"+piece_size+"-"+date+".dat",'w')
order_h = open("order-"+piece_size+"-"+date+".dat",'w')
class_h = open("class-"+piece_size+"-"+date+".dat",'w')
phylum_h = open("phylum-"+piece_size+"-"+date+".dat",'w')

token_split_tab = re.compile(r'\t')
for line in file_h:
    tokens = token_split_tab.split(line)

    print 'Fetching: %s' % str(tokens)
    cur_spec_tax = fetchTaxonomy(tokens[1], cursor)
    high_spec_tax = fetchTaxonomy(tokens[2], cursor)

    if cur_spec_tax == None or high_spec_tax == None: continue

    if 'genus' in high_spec_tax.keys() and 'genus' in cur_spec_tax.keys():
        genus_h.write(tokens[0]+'\t'+cur_spec_tax['genus']+'\t'+
                      high_spec_tax['genus']+'\t'+tokens[3])
    if 'family' in high_spec_tax.keys() and 'family' in cur_spec_tax.keys():
        family_h.write(tokens[0]+'\t'+cur_spec_tax['family']+'\t'+
                       high_spec_tax['family']+'\t'+tokens[3])
    if 'order' in high_spec_tax.keys() and 'order' in cur_spec_tax.keys():
        order_h.write(tokens[0]+'\t'+cur_spec_tax['order']+'\t'+
                      high_spec_tax['order']+'\t'+tokens[3])
    if 'class' in high_spec_tax.keys() and 'class' in cur_spec_tax.keys():
        class_h.write(tokens[0]+'\t'+cur_spec_tax['class']+'\t'+
                      high_spec_tax['class']+'\t'+tokens[3])
    if 'phylum' in high_spec_tax.keys() and 'phylum' in cur_spec_tax.keys():
        phylum_h.write(tokens[0]+'\t'+cur_spec_tax['phylum']+'\t'+
                       high_spec_tax['phylum']+'\t'+tokens[3])

file_h.close
genus_h.close
family_h.close
order_h.close
class_h.close
phylum_h.close

db_h.close
print "DONE DOIN' WORK SON!"
