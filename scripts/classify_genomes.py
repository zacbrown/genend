#!/usr/bin/env python

import MySQLdb
import yaml
import sys

if len(sys.argv) < 4:
    print "ERROR while doin' work son!"
    sys.exit(-1)

# load the config
config_path = sys.argv[3]
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

        query = "SELECT name FROM taxon_name WHERE "
        query += "name_class='scientific name' AND "
        query += "taxon_id='%s'" % cur_taxon_id

        cursor.execute(query)

        ret_dict[rank] = cursor.fetchone()[0]

        cur_taxon_id = next_taxon_id

    return ret_dict




# load file
file_h = open(sys.argv[1], "rU")

db_h = MySQLdb.connect(host=db_vals["server"], user=db_vals["username"],
                       passwd=db_vals["password"], db=db_vals["database"])

match_vals = {"kingdom":0, "phylum":0, "class":0, "order":0, "family":0, "genus":0}
total_vals = {"kingdom":0, "phylum":0, "class":0, "order":0, "family":0, "genus":0}

cursor = db_h.cursor()
count_measured = 0

for line in file_h:
    (cur_spec, sep, high_spec) = line.partition("\t")

    print 'Fetching: %s | %s' % (cur_spec, high_spec)
    cur_spec_tax = fetchTaxonomy(cur_spec, cursor)
    high_spec_tax = fetchTaxonomy(high_spec, cursor)

    if cur_spec_tax == None or high_spec_tax == None: continue

    for key in match_vals:
        if key in cur_spec_tax and key in high_spec_tax:
            total_vals[key] += 1
            if cur_spec_tax[key] == high_spec_tax[key]: match_vals[key] += 1

file_h.close
output_h = open(sys.argv[2], 'w')

print match_vals

for key, val in match_vals.iteritems():
    output_h.write(key + "\t" + str(val) + "\t" + str(total_vals[key]) + "\n")

output_h.close
db_h.close
print "DONE DOIN' WORK SON!"
