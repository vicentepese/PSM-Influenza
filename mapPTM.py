from Bio import Entrez
from Bio import SeqIO
import numpy as np 
import sys 
import os 
import json 
import csv 
from collections import defaultdict, OrderedDict
import re 
from markdown2 import Markdown

def importData(options):

	# Import data  as list of lists
	data = list() 
	with open(options['files']['mergedFiltMASS'], 'r') as inFile:
		reader = csv.reader(inFile)
		for row in reader:
			data.append(row)
	
	return data

def reference_retreive(proteinID):

	# Retrieve protein and get dictionnary of each peptide position
	Entrez.email = "ambati@stanford.edu"
	handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=proteinID)
	seq_record = SeqIO.read(handle, "fasta")

	# For each reference
	out_ref={}
	for i, j in enumerate(str(seq_record.seq)):
		out_ref[i+1] =j
	return out_ref

def map_PTMs(options, data, refProt):

    # Initialize 
    PTM_map = defaultdict()
    for pos in list(refProt.keys()):
        PTM_map[pos] = defaultdict(lambda: defaultdict(int))

    # For each fragment
    for seq in data:

        # Get initial position 
        init_pos = int(seq[2])
        AAseq = seq[1][2:-2]

        # If initial position in range and there is a PTM
        if init_pos >= options['pos_range'][0] and init_pos <= options['pos_range'][1] and '[' in seq[1]:

            # Get PTM type and location 
            PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
            PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

            # For each PTM
            idx_cumm = 0
            for instance, idx in zip(PTM_instances, PTM_idx):

                # Find position 
                ptm_pos = init_pos + idx.start() - 1 - idx_cumm
                PTM_map[ptm_pos][instance][seq[3]] += 1
                idx_cumm += len(instance) + 2
    
    return PTM_map

def map2HTML(PTM_map):

    PTM_mark = list()
    # For each position, PTM, and vaccine 
    for pos in list(PTM_map.keys()):
        for ptm in list(PTM_map[pos].keys()):
                PTM_mark.append([str(pos) + ': ' + str(ptm) + '(' + 'PAN:' + str(PTM_map[pos][ptm]['PAN']) + ' ARP:' + str(PTM_map[pos][ptm]['ARP']) + ')' + '\n'])

    # Convert to HTML
    PTM_HTML = list()
    markdowner = Markdown()
    for mark_sent in PTM_mark:
        PTM_HTML.append(markdowner.convert(mark_sent))
    
    # Write 
    with open('test.html','r') as outFile:
        outFile.write(sent for sent in PTM_HTML)

                
def main():

    # Read options 
    with open('options.json', 'r') as inFile:
        options = json.load(inFile)

    # Import data
    data = importData(options)

    # Import protein of reference
    refProt = reference_retreive(options['refProt'])

    # Map PTMs 
    PTM_map = map_PTMs(options, data, refProt)
    
    # Convert to HTML and store
    map2HTML(PTM_map)


if __name__ == "__main__":
    main() 
