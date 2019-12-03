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
    vaccSample = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    # For each fragment
    for seq in data:

        # Get initial position and sequence of AA
        init_pos = int(seq[2])
        AAseq = seq[1][2:-2]

        # Count 
        AAnonPTM = re.sub('\[.+?\]','',seq[1]).split('.')[1]
        end_pos = init_pos + len(AAnonPTM)
        for i in range(int(seq[2]), int(seq[2]) + len(AAnonPTM)):
            vaccSample[str(i)][AAnonPTM[i-init_pos]][seq[3]] += 1

        # If initial position in range and there is a PTM
        if not(end_pos < options['pos_range'][0]) and not(init_pos > options['pos_range'][1]) and  '[' in seq[1]:

            PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
            PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

            # For each PTM
            idx_cumm = 0
            for instance, idx in zip(PTM_instances, PTM_idx):

                # Find position 
                ptm_pos = init_pos + idx.start() - 1 - idx_cumm
                PTM_map[ptm_pos][instance][seq[3]] += 1
                idx_cumm += len(instance) + 2
                
    return PTM_map, vaccSample

def map2HTML(PTM_map, refProt, vaccSample, options):

    # For each position, PTM, and vaccine 
    PTM_HTML = list()
    markdowner = Markdown()
    for pos in range(options['pos_range'][0], options['pos_range'][1]+1):
        if len(list(PTM_map[pos].keys())) >= 1:
            PTM_mark = str(refProt[pos]) + ': ' 
            for ptm in list(PTM_map[pos].keys()):
                PTM_mark = PTM_mark + '__' + str(ptm) + '__' + '(ARP:{:.2%}' ' PAN:{:.2%}\) &nbsp;'.format(PTM_map[pos][ptm]['ARP']/vaccSample[str(pos)][refProt[pos]]['ARP'], PTM_map[pos][ptm]['PAN']/vaccSample[str(pos)][refProt[pos]]['PAN'])
            PTM_mark = PTM_mark + ' \n'
            PTM_HTML.append(markdowner.convert(PTM_mark))
        else:
            PTM_mark = str(refProt[pos]) + ' \n'
            PTM_HTML.append(markdowner.convert(PTM_mark))

    # Write 
    with open(options['files']['mapPTM.html'],'w') as outFile:
        # Header defining html style 
        outFile.write('<style>' + '\n' + 'p { \n' + 'line-height:0.01; \n' + 'font-family: "Courier New", Courier, monospace; \n' + '} \n' + '\n' +"</style>" +'\n')
        outFile.write(markdowner.convert(str(options['pos_range'][0]) + '\n'))
        outFile.writelines(PTM_HTML)
        outFile.write(str(options['pos_range'][1]))

                
def main():

    # Read options 
    with open('options.json', 'r') as inFile:
        options = json.load(inFile)

    # Import data
    data = importData(options)

    # Import protein of reference
    refProt = reference_retreive(options['refProt'])

    # Map PTMs 
    PTM_map, vaccSample = map_PTMs(options, data, refProt)
    
    # Convert to HTML and store
    map2HTML(PTM_map, refProt, vaccSample, options)


if __name__ == "__main__":
    main() 
