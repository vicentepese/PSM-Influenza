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
from scipy import stats 

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

def div0(n, d):
	return n / d if d and n else 0

def map_PTMs(options, data, refProt):

    # Initialize 
    PTM_map = defaultdict()
    for pos in refProt:
        PTM_map[pos] = defaultdict(lambda: defaultdict(int))
    vaccSample = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    # For each fragment
    for seq in data:

        # Get initial position and sequence of AA
        AAseq = seq[1][2:-2]
        AAnonPTM = re.sub('\[.+?\]','',AAseq)
        init_pos = int(seq[2])
        end_pos = init_pos + len(AAnonPTM)

        # If initial position in range and there is a PTM
        if not(end_pos < options['pos_range'][0]) and not(init_pos > options['pos_range'][1]):

            if '[' in seq[1]:
                PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
                PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

                # For each PTM, find position and append according to type 
                idx_cumm = 0
                for instance, idx in zip(PTM_instances, PTM_idx):
                    ptm_pos = init_pos + idx.start() - 1 - idx_cumm
                    PTM_map[ptm_pos][instance][seq[3]] += 1
                    idx_cumm += len(instance) + 2
            
            # Count 
            for i in range(int(seq[2]), int(seq[2]) + len(AAnonPTM)):
                vaccSample[i][AAnonPTM[i-init_pos]][seq[3]] += 1

                
    return PTM_map, vaccSample

def statisticalTest(options, PTM_map, vaccSample, refProt):

    # Initialize 
    PTM_stats = defaultdict(lambda: defaultdict(lambda : defaultdict(int)))

    # For each position
    for pos in range(options['pos_range'][0], options['pos_range'][1]+1):

        # IF there is a ptm for both vaccines
        if len(list(PTM_map[pos].keys())) >=1:
            for ptm in list(PTM_map[pos].keys()):
                if PTM_map[pos][ptm]['PAN'] and PTM_map[pos][ptm]['ARP']:

                    # Create array 
                    ptm_positive = [PTM_map[pos][ptm]['ARP'], PTM_map[pos][ptm]['PAN']]
                    ptm_negative = [vaccSample[pos][refProt[pos]]['ARP'] - PTM_map[pos][ptm]['ARP'], \
                         vaccSample[pos][refProt[pos]]['PAN'] - PTM_map[pos][ptm]['PAN']]
                    
                    # Fisher test and append to output
                    oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
                    PTM_stats[pos][ptm]['pvalue'] =  pvalue
                    PTM_stats[pos][ptm]['oddsratio'] = oddsratio


    return PTM_stats

def map2HTML(PTM_map, refProt, vaccSample, options, PTM_stats):

    # For each position, PTM, and vaccine 
    PTM_HTML = list()
    markdowner = Markdown()
    for pos in range(options['pos_range'][0], options['pos_range'][1]+1):
        if len(list(PTM_map[pos].keys())) >= 1:
            PTM_mark = str(refProt[pos]) + ': (ARP: {0}, PAN:{1}) // '.format(vaccSample[pos][refProt[pos]]['ARP'], vaccSample[pos][refProt[pos]]['PAN']) 
            for ptm in list(PTM_map[pos].keys()):

                if not PTM_stats[pos][ptm]:
                    PTM_mark = PTM_mark + '__' + str(ptm) + '__' + \
                        '(ARP:{:.2%}' ' PAN:{:.2%}) &nbsp;'.format(div0(PTM_map[pos][ptm]['ARP'],vaccSample[pos][refProt[pos]]['ARP']),\
                            div0(PTM_map[pos][ptm]['PAN'],vaccSample[pos][refProt[pos]]['PAN']))
                elif PTM_stats[pos][ptm]['pvalue'] > 0.05:
                    PTM_mark = PTM_mark + '__' + str(ptm) + '__' + \
                        '(ARP:{:.2%}' ' PAN:{:.2%}, p = {:.2}) &nbsp;'.format(div0(PTM_map[pos][ptm]['ARP'],vaccSample[pos][refProt[pos]]['ARP']),\
                            div0(PTM_map[pos][ptm]['PAN'],vaccSample[pos][refProt[pos]]['PAN']), PTM_stats[pos][ptm]['pvalue']) 
                else:
                    PTM_mark = PTM_mark + '__' + str(ptm) + '__' + \
                        '(ARP:{:.2%}' ' PAN:{:.2%},  <span style=\"color: red;\"> p = {:.2}</span>) &nbsp;'.format(div0(PTM_map[pos][ptm]['ARP'],vaccSample[pos][refProt[pos]]['ARP']),\
                            div0(PTM_map[pos][ptm]['PAN'],vaccSample[pos][refProt[pos]]['PAN']), PTM_stats[pos][ptm]['pvalue']) 
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

    # Statistical test 
    PTM_stats = statisticalTest(options, PTM_map, vaccSample, refProt)
    
    # Convert to HTML and store
    map2HTML(PTM_map, refProt, vaccSample, options, PTM_stats)


if __name__ == "__main__":
    main() 
