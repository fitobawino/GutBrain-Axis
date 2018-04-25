# Uses the GCRh38 version of RefSeq proteins to create a table of accession to geneID file.

import csv

file = '/.../ftp.ncbi.nlm.nih.gov/interim_GRCh38.p10.RefSeqProts.gpff'
save_path = '/.../acc2geneID.csv'

dat = dict()
with open(file) as lines:
	ver = False
	for line in lines:
		if line.split() == []: continue
		if line.split()[0] == 'VERSION':
			acc = line.split()[-1]
			ver = True
		if ver and len(line.split()[0])>17:
			if line.split()[0][0:17] == '/db_xref="GeneID:':
				dat[acc] = line.split()[0][17:-1]
				ver = False

with open(save_path,'wb') as file:
	writer = csv.writer(file)
	for key,value in dat.items():
		writer.writerow([key, value])