from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW
from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import semantic_similarity, TermCounts
import os, datetime, csv, sys
from os.path import join
import pandas as pd
import numpy as np 

def HMGvsHum(seq_file, refseq_db, blast_res, save_path):
# Psi-Blast metagenome vs human
# 'seq_file' is the path to fasta-formated protein sequences
# 'refseq_db' is the path to the database to BLAST to.
# 'blast_Res' is the save path
# 'save_path' is the path where the psiblast results will be saved.

	os.system("touch " + join(blast_res, save_path))
	result_handle = os.system("psiblast" + 
    	" -num_threads 8" +
    	" -num_iterations 10" +
    	" -evalue 0.05" +
    	" -db " + refseq_db +
    	" -outfmt 7" +
    	" -query " + seq_file +
    	" -out " + join(blast_res, save_path))


def Blast2genes(blast_file, save_path, path_acc2id, taxonomy):
# Program that parses the psi-blast result in the path 'blast_file', 
# retrieve genes from the association file 'path_acc2id' created with 
# Acc2geneIDdb.py and saves the list to 'save_path'

    df        = pd.DataFrame({'H.prot':[], 'Tax':[]})
    prot_ids  = list()
    tax       = taxonomy
    it10      = False
    count     = 1
    with open(blast_file) as lines:
        reading = False
        for line in lines:
            words = line.split()
            if words != []: 
                if line.startswith('# Iteration: '): #Identify iter=10
                    if line[13:] == '10\n': 
                        it10 = True
                    else: it10 = False
                # Erase prot_temp var
                if words == ['#', 'PSIBLAST','2.6.0+']:
                    if it10:
                        for prot in prot_temp: prot_ids.append(prot)
                    prot_temp = list()
                    reading = False  
                # begin to save proteins list
                if reading: prot_temp.append(line.split()[1]) #Write proteinID to var            
                if words[1] != '0':
                    if words[-2:] == ['hits', 'found']: reading = True
            if line == '\n':
                for prot in prot_temp: prot_ids.append(prot)
                prot_ids = list(set(prot_ids))
                prot_temp = list()
                it10 = False
    prot_ids = list(set(prot_ids)) #Remove duplicates
    #Create dataframe
    df = pd.DataFrame.from_items([('H.prot', prot_ids)])
# Change prot ref to geneID
    geneID = list()
    with open(path_acc2id, 'r') as csv_file:
        reader = csv.reader(csv_file)
        acc2id = dict(reader)
    for prot_id in prot_ids:
        if prot_id in acc2id: 
            geneID.append(acc2id[prot_id])
        else: 
            geneID.append('')
    df['geneIDs'] = pd.Series(geneID)
#save df to file
    df.to_csv(join(save_path, 'df_' + tax.replace(' ', '_') + '.csv'))
#Write the genes list to a file
    geneID = list(set(geneID)) #Remove duplicates
    if geneID[0] == '': geneID.remove('') #Remove empty values
    geneID_route = join(save_path, 'geneIDs_' + tax.replace(' ', '_') + '.txt')
    with open(geneID_route, 'w') as myfile:
        for line in zip(geneID):
            myfile.write("{0}\n".format(*line))


def genes2GO(study_path, pop_path, association_path, go_db_path, save_path, filename, overonly):
# function that uses the geneIDs and performs GO enrichment.
# 'study_path' path to a file of geneIDs in one column
# 'pop_path' is the path to a file with all of the geneIDs in the population (same format).
# 'associaton_path' is the path to a file that with two columns: <geneID> <GO1;GO2;...GOn> (created by gene2GOdb.py)
# 'save_path' is the path to where the enrichment will be saved

#Load study genes
    os.system('find_enrichment.py '+ 
    '--pval=0.05 --alpha=0.05 ' +
    '--obo ' + go_db_path + ' ' +
    study_path + ' ' +
    pop_path + ' ' +
    association_path + ' ' + 
    #'--pvalcalc fisher '
    '--outfile ' + join(save_path, "aux_" + filename))
#Bonferroni-threshold
    with open(join(save_path, "aux_" + filename), 'r+') as aux_file, open(join(save_path, filename), 'w') as cvs_file:
        reader = csv.reader(aux_file,  delimiter='\t')
        next(reader, None) 
        writer = csv.writer(cvs_file,  delimiter='\t')
        for row in reader:
             if (float(row[9]) <= 0.05):
             	if ((overonly) & (row[2] == 'e')) | (not overonly):
             		writer.writerow(row)


def filter_GO_by_slim(en_path, slim_path, save_path):
# This function removes the ontologies from the enrichment file that do not belong to the 
# GO_slim

    go_list = list()
    with open(slim_path,'r') as slim_file:
        reader = csv.reader(slim_file, delimiter = '\t')
        for line in reader:
            go_list.append(line[4])
    go_inter = list()
    with open(en_path, 'r') as en_file:
        reader = csv.reader(en_file, delimiter = '\t')
        for line in reader: 
            if line[0] in go_list:
                if float(line[9]) <= 0.05:
                    go_inter.append(line)
# if its empty do not write
    if len(go_inter) > 0:
        with open(save_path,'w') as save_file:
            wr = csv.writer(save_file, delimiter = '\t')
            for row in go_inter:
                wr.writerow(row)


def find_GO_intersection(path1, path2, gene2des_path, starts_with, save_path):
#path1 and path2 are paths to files starting with 'starts_with' that are tsv with GO 
#data at the first row. This function compares all GOs in path1 with all GOs in path2

    path1_files = [files for files in os.listdir(path1) if files.startswith(starts_with)]
    path2_files = [files for files in os.listdir(path2) if files.startswith(starts_with)]
    prefix = 'inter_'
    #Creates a dictionary with GO -> description.
    go2description  = dict()
    csv_file = open(gene2des_path, 'r')
    reader = csv.reader(csv_file)
    for row in reader: go2description[row[2]] = row[5]
    for file1 in path1_files:
        go_list1 = list()
        file = open(join(path1, file1),'r')
        reader = csv.reader(file, delimiter = '\t')
        for line in reader: go_list1.append(line[0])
        for file2 in path2_files:
            go_list2 = list()
            file = open(join(path2, file2), 'r')
            reader = csv.reader(file, delimiter = '\t')
            for line in reader: 
            	go_list2.append(line[0])
            inter = set(go_list1) & set(go_list2)
            inter = list(inter)
            #String with descriptions separated by ';'
            descr = str()
            inter = [GOs for GOs in inter if GOs in go2description]
            if len(inter) > 0:
                res = str()
                save_path2 = join(save_path, 'results_all.txt')
                for elem in inter: res = res + ',' + elem 
                with open(save_path2, 'a') as save_file:
                    wr = csv.writer(save_file, delimiter = '\t')
                    wr.writerow([file1[len(starts_with):-4], file2[len(starts_with):-4], res[1:], descr])
                save_path3 = join(save_path, prefix) + file1[len(starts_with):-4] + '-WITH-' + file2[len(starts_with):-4] + '.txt'
                with open(save_path3, 'a') as save_file:
                    for line in zip(inter):
                        save_file.write("{0}\t".format(*line))
                        save_file.write(go2description[line[0]] + "\n")


def merge_by_genus(files_path, save_path):
# Remove redundancies from several files and creates a single fasta file

    pep_files = os.listdir(files_path)
    nr_seqs = []
    for pep_file in pep_files:
        print("  " + pep_file)
        file_seqs = list(SeqIO.parse(join(files_path, pep_file), 'fasta'))
        if len(nr_seqs) == 0:
            nr_seqs = file_seqs
            continue
        for file_seq in file_seqs: #Seqs in file
            for nr_seq in nr_seqs: #Non redundant seqs
                if nr_seq.seq == file_seq.seq:
                    novel = False
                    break
                else:
                    novel = True
            if novel: nr_seqs.append(file_seq)
    SeqIO.write(nr_seqs, save_path, "fasta")


def GO_distance(file_path, save_path, go_db):
# Finds the distance matrix of the intersection GOs and saves it.

    gos = list()
    prefix = 'dist_'
    with open(file_path, 'r') as csv_file:
        reader = csv.reader(csv_file,  delimiter='\t')
        for row in reader: gos.append(row[0])
    dist_mat = np.matrix([[0.0000] * len(gos)]*len(gos))
    for num in range(0,len(gos)):
        for num2 in range(num + 1,len(gos)):
            try:
                dist_mat[num, num2] = semantic_similarity(gos[num], gos[num2], go_db)
                #if np.isnan(dist_mat[num, num2]): dist_mat[num, num2] = 0
                dist_mat[num2, num] = dist_mat[num, num2]
            except:
                print(num, num2)
    np.savetxt(save_path, dist_mat)
    return dist_mat
