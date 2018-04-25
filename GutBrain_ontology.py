'''
Script used for the following work:
AUTHOR=Flores Saiffe Farías Adolfo, Mendizabal Adriana P., Morales J. Alejandro
TITLE=An Ontology Systems Approach on Human Brain Expression and Metaproteomics
JOURNAL=Frontiers in Microbiology     
VOLUME=9      
YEAR=2018
PAGES=406   
URL=https://www.frontiersin.org/article/10.3389/fmicb.2018.00406     
DOI=10.3389/fmicb.2018.00406    
ISSN=1664-302X   

Created by Adolfo Flores Saiffe Farías
Departamento de Ciencias Computacionales
Centro Universitario de Ciencias Exactas e Ingenierías
Universidad de Guadalajara
adolfo.flores.saiffe@gmail.com
saiffe.adolfo@alumnos.udg.mx

Requirements:
python version >= 3.5
psiblast

Packages:
 - Biopython
 - goatools
 - pandas

'''

# Import packages
import sys, os
utils_dir = '/.../FW_test/Scripts/'
sys.path.append(utils_dir)
import GutBrain_ontology_utils as GBOu
from os.path import join
from itertools import compress
from goatools import obo_parser

#################
# Set variables #
#################
#1. Merge .pep.fsa files per genre
working_dir   = '/../FW_test'
databases_dir = join(working_dir, 'Databases')
results_dir   = join(working_dir, 'Results')
meta_route 	  = join(databases_dir, 'HMP')

#2. PSI-Blast and 3. Parsing psi-blast results into list of genes
refseq_path   = "/.../Downloads/ncbi-blast-2.6.0+/bin/refseq_prot"
blast_subdir  = 'Blast'

#4. Enrichment of microbiota-gene and filtering by Bonferroni < 0.05
enrich_dir    = 'GO-enrichment'
path_acc2id   = join(databases_dir, 'acc2geneID.csv') # File created by GutBrain_Acc2geneIDdb.py
pop_path      = join(databases_dir, 'pop.txt') # File created by GutBrain_gene2GOdb.py
gene2go_path  = join(databases_dir, 'gene2go.txt') # File created by GutBrain_gene2GOdb.py
go_db_path    = join(databases_dir, "GeneOntology/go-basic.obo")

#5. Enrichment of brain-gene expression and filtering by Bonferroni < 0.05 
# Path to genes differentially expressed genes (output of Exp_data_analysis.r)
exp_genes_dir = join(results_dir, 'Expression')

# 6. Find the intersection between the two enrichments,
inter_save_path= join(results_dir, "GOs-intersection")
go_db_path     = join(databases_dir, "GeneOntology/go.obo")
dist_save_path = join(results_dir, "Distance")

# 7. Create table S1: Taxon | Metaproteome | Proteins found by psiblast | no. of GOs in enrichment
T1_save_path = join(results_dir, "Tables/TableS1.csv")

# 8. Create Table S2: Substructure | no. genes with diff.exp. | others | no of GOs in enrichment
T3_save_path = join(results_dir, "Tables/TableS3.csv")

os.mkdir(working_dir)
os.mkdir(databases_dir)
os.mkdir(results_dir)
os.mkdir(meta_route)
os.mkdir(join(results_dir, blast_subdir))
os.mkdir(join(results_dir, enrich_dir))
os.mkdir(exp_genes_dir)
os.mkdir(dist_save_path)
os.mkdir(inter_save_path)
os.mkdir(join(results_dir, 'Tables'))

##############
# Processing #
##############

# 1. Merge .pep.fsa files by genre
genres  = os.listdir(meta_route)
for genre in genres:
	print("Merging files from: " + genre)
	GBOu.merge_by_genus(
		files_path = join(meta_route, genre),
		save_path = join(meta_route, genre + '.pep.fsa'))

# 2. PSI-Blast study: Loop to psi-blast all files in seq_files.
query_files = [x for x in os.listdir(meta_route) if x.endswith('.fsa')]
for query_file in query_files:
	print ("Query file: " + query_file)
	GBOu.HMGvsHum(
		seq_file   = join(meta_route, query_file),
		refseq_db  = refseq_path,
		blast_res  = join(results_dir, blast_subdir),
		save_path  = query_file[:-8])

# 3. Parsing psi-blast results into list of genes ('df_<taxon>.csv', 'geneIDs_<taxon>.txt').
blast_files = os.listdir(join(results_dir, blast_subdir))
for file in blast_files:
	print ("Processing: " + file)
	GBOu.Blast2genes(
        blast_file  = join(results_dir, blast_subdir, file),
        taxonomy    = file,
        save_path   = join(results_dir, enrich_dir),
        path_acc2id = path_acc2id)

# 4. Enrichment of microbiota-gene and filtering by Bonferroni < 0.05 ('GOen_genus_<taxon>')
genes_path = join(results_dir, enrich_dir)
gene_files = [filename for filename in os.listdir(genes_path) if filename.startswith("geneIDs_")]
for gene_file in gene_files:
	print ("\nProcessing: " + gene_file)
	GBOu.genes2GO(
		study_path       = join(genes_path, gene_file),
		pop_path         = pop_path,
		association_path = gene2go_path,
		go_db_path 		 = go_db_path,
		save_path        = genes_path,
		filename         = 'GOen_genus_' + gene_file[8:],
		overonly         = True)

# 5. Enrichment of brain-gene expression and filtering by Bonferroni < 0.05 
# Genes with no annotation will drop a warning message.
gene_files = [filename for filename in os.listdir(exp_genes_dir) if filename.startswith("Genes_logFC")] #cuales usar?
for gene_file in gene_files:
	print ("\nProcessing: " + gene_file)
	GBOu.genes2GO(
		study_path       = join(exp_genes_dir, gene_file), 
		pop_path         = pop_path, 
		association_path = gene2go_path, 
		go_db_path 		 = go_db_path,
		save_path        = exp_genes_dir,
		filename         = 'GOen_genus_' + gene_file[8:],
		overonly         = False) #revisar si funciona bien

# 6. Find the intersection between the two enrichments, calculate the distance matrix and saves it
genes_path = join(results_dir, enrich_dir)
prefix = 'GOen_genus_'
GBOu.find_GO_intersection(
	path1        = genes_path, 
	path2        = exp_genes_dir,
	gene2des_path = join(databases_dir, 'GeneOntology', 'gene2go_human.csv'),
	starts_with  = prefix,
	save_path    = inter_save_path)

go_db = obo_parser.GODag(go_db_path)
files = [files for files in os.listdir(inter_save_path) if files.startswith('inter_')]
for file in files:
	GBOu.GO_distance(
		file_path = join(inter_save_path, file),
		save_path = join(dist_save_path, "dist_" + file),
		go_db     = go_db)

# 7. Create table S1: Taxon | Metaproteome | Proteins found by psiblast | no. of GOs in enrichment
from Bio import SeqIO
import pandas as pd
import json
Bprots = list()
Hprots = list()
GOs    = list()
metaproteome_files = [x for x in os.listdir(meta_route) if x.endswith('.fsa')]
for taxon in metaproteome_files:
	Bprots.append(len(list(SeqIO.parse(join(meta_route, taxon), "fasta"))))
	with open(join(results_dir, enrich_dir, 'geneIDs_' + taxon[:-8] + '.txt'), 'r') as myfile:
		for Hprot, l in enumerate(myfile):
			pass
		Hprots.append(Hprot + 1)

	with open(join(results_dir, enrich_dir, 'GOen_genus_' + taxon[:-8] + '.txt'), 'r') as myfile:
		for GO, l in enumerate(myfile):
			pass
		GOs.append(GO + 1)

todos = {'Taxon': metaproteome_files,
		 'Metaproteome': Bprots,
		 'Prots_found_by_blast': Hprots,
		 'GOs_in_enrichment': GOs}
df = pd.DataFrame.from_dict(todos)
df.to_csv(path_or_buf = T1_save_path, sep = '\t')

# 8. Create Table S3: Substructure | no. genes with diff.exp. | others | no of GOs in enrichment
Hprots = list()
GOs    = list()
subs_files = [x for x in os.listdir(exp_genes_dir) if x.startswith('Genes_logFC')]
for sub in subs_files:
	with open(join(exp_genes_dir, sub), 'r') as myfile:
		for Hprot, l in enumerate(myfile):
			pass
		Hprots.append(Hprot + 1)
	with open(join(exp_genes_dir, 'GOen_genus_gFC.' + sub[12:]), 'r') as myfile:
		for GO, l in enumerate(myfile):
			pass
		GOs.append(GO + 1)
todos = {'Substructures': subs_files,
		 'Diff._exp._genes': Hprots,
		 'GOs_in_enrichment': GOs}
df = pd.DataFrame.from_dict(todos)
df.to_csv(path_or_buf = T3_save_path, sep = '\t')

# 9. Find the distances between all GOs
go_db = obo_parser.GODag(go_db_path)
GBOu.GO_distance(
	file_path = join(databases_dir, "GeneOntology/GO2label_handcurated.csv"),
	save_path = join(results_dir, "Distance/All_GOs.csv"),
	go_db     = go_db)