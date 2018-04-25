# Script that imports gene2go_human.csv and creates an association file with the following format:
# <geneID1> <GO1;GO2>
# <geneID2> <GO3;GO4;GO5>
# ...
# <geneIDn> <GOm>
# Also creates a list of human genes by geneid and saves it at pop_save_path
# gene2go_human.csv was mannually created from the file: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

#Load libraries
import csv
import numpy as np

#Variables
#Path to gene2go_human.csv
import_path   = '/.../gene2go_human.csv'
#Save path to association file
g2g_save_path = '/.../gene2go.txt'
#Save path to human gene list
pop_save_path = '/.../pop.txt'
gene2go  = dict()

#Import file and merge GOs from the same geneID.
with open(import_path, 'rb') as csv_file:
    reader = csv.reader(csv_file)
    for row in reader:
    	if gene2go.has_key(row[1]):
    		gene2go[row[1]] = gene2go[row[1]] + ';' + row[2]
    	else: gene2go[row[1]] = row[2]

#Save var to file
keys = gene2go.keys()
with open(g2g_save_path, 'w') as myfile:
    for key in keys:
        myfile.write(key + '\t' + gene2go[key] + '\n')

#Create a file with all gene IDs
with open(pop_save_path, 'w') as myfile:
    for key in keys:
        myfile.write(key + '\n')