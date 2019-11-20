from Bio import Entrez
from Bio import SeqIO
import json
from subprocess import call
import csv
from collections import defaultdict
import numpy as np 
import re


## import filenames as a list
f_names=['ar_chy_list','ar_try_list', 'pan_chy_list', 'pan_try_list']### def subroutines

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

def filterData(options):
    
    # Filt in bash 
    call(options['filterDataBash'], shell=True)

def importAFLSA(options):

    # Import data  as list of lists
    AFLSA = list() 
    with open(options['files']['filt_AFLSA'], 'r') as inFile:
        reader = csv.reader(inFile)
        for row in reader:
            AFLSA.append(row)
    
    return AFLSA

def mapSeqs(options, AFLSA, refProt):

    # Get proteins 
    proteins = list(np.unique(np.asarray([seq[0] for seq in AFLSA])))

    # For each sequence 
    pos_range = options['pos_range']
    seqMap = defaultdict(list)
    for seq in AFLSA:

        # Get initial position
        init_pos = int(seq[-1])

        # Remove PTM and compute ending positon, else no PTM and compute ending position
        if '[' in seq[1]:
            
            # Get sequence and ending position
            AAseq = re.sub('\[.+?\]','',seq[1]).split('.')[1]
            end_pos = init_pos  + len(AAseq)-1
        else:
            # Get sequence
            AAseq = seq[1].split('.')[1]
            end_pos = init_pos + len(AAseq)-1

        # If initial position in the range, else check if any position overlap
        if init_pos >= pos_range[0] and init_pos <= pos_range[1]:

            # Map each position 
            for pos, AA in zip(range(init_pos, (end_pos + 1)), AAseq):
                if pos == 250 and AA == "L":
                    print(AAseq)
                seqMap[str(pos)].append(AA)
        
        elif end_pos >= pos_range[0] and end_pos <= pos_range[1]:

            # Map each position 
            for pos, AA in zip(range(init_pos, (end_pos + 1)), AAseq):
                # 'HA-1189-G (397) ((D)) (((N)))'
                if pos == 250 and AA == "Q":
                    print(seq[0])
                seqMap[str(pos)].append(AA)
            
    return seqMap

def countSequences(options, AFLSA):

    # Get proteins 
    proteins = list(np.unique(np.asarray([seq[0] for seq in AFLSA])))

    # For each sequence 
    pos_range = options['pos_range']
    seqMap = defaultdict(list)

    # TODO: Count sequences

    pass



def main():

    # Read options
    with open('options.json','r') as inFile:
        options = json.load(inFile)
    
    # Filter data 
    filterData(options)

    # Get reference protein 
    refProt = reference_retreive("ACQ55359.1")

    # Import AFLSA
    AFLSA = importAFLSA(options)

    # Map sequences 
    seqMap = mapSeqs(options, AFLSA, refProt)

    for pos in list(seqMap.keys()):
        seqMap[pos] = list(np.unique(np.asarray(seqMap[pos])))

    # Count sequences
    countSequences(options, AFLSA)

if __name__ == "__main__":
    main()