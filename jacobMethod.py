import numpy as np 
import sys 
import os 
import json 
import csv 
import re
import random 
import subprocess
from markdown2 import Markdown
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict, OrderedDict
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

def importBindData(options):

	# Import binders 
	binders = defaultdict(int)
	with open(options['files']['DQ0602Binders'],'r') as inFile:
		next(inFile)
		reader = csv.reader(inFile)
		for row in reader:
			binders[row[2]] = float(row[3])
	
	# Import binding cores 
	bindingCores = defaultdict(str)
	with open(options['files']['DQ0602BindingCores'], 'r') as inFile:
		next(inFile)
		reader = csv.reader(inFile, delimiter = '\t')
		for row in reader:
			bindingCores[row[4]] = row[6]

	return binders, bindingCores


def getRandomColor(options, PTM_count):

	# Create color dictionnary	
	r = lambda: random.randint(75,200)
	color = {PTM: ['<span style=\"background: '+'#%02X%02X%02X; font-weight: bold' % (r(),r(),r()) + '\">',"</span>"] for PTM in list(PTM_count.keys())}
	color['red'] = ['<span style=\"background: red; font-weight: bold; \">', '</span>']
	color['orange'] = ['<span style=\"background: orange; font-weight: bold; \">', '</span>']
	color['bindingCore'] = ['<span style=\"background: green; font-weight: bold; \">', '</span>']

	return color

def div0(n, d):
	return n / d if d and n else 0


def countPositions(data):

	# Position count 
	posCount = defaultdict(lambda: defaultdict(int))

	# Count number of occurences of an AA for each position 
	for seq in data:
		
		# Clear data and count AA 
		AAseq = re.sub('\[.+?\]','', seq[1][2:-2])
		for i in range(0,len(AAseq)):
			posCount[int(seq[2])+i][AAseq[i]] += 1 

	return posCount

def statisticalTest(options, PTM_seq, seqCount, seq):

	# Initialize 
	PTM_stats = defaultdict(lambda: defaultdict(lambda : defaultdict(int)))
	sig_pval = list()

	# For each ptm
	for pos in list(PTM_seq.keys()):
		for ptm in list(PTM_seq[pos].keys()):
			if PTM_seq[pos][ptm]['PAN'] and PTM_seq[pos][ptm]['ARP']:

				# Create array
				ptm_positive = [PTM_seq[pos][ptm]['ARP'], PTM_seq[pos][ptm]['ARP']]
				ptm_negative = [seqCount[seq]['ARP'] - PTM_seq[pos][ptm]['ARP'], \
					seqCount[seq]['PAN'] - PTM_seq[pos][ptm]['PAN']]
				
				# Fisher test, append to output 
				oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
				PTM_stats[pos][ptm]['pvalue'] = pvalue
				PTM_stats[pos][ptm]['oddsratio'] = oddsratio

				# Append significant pvalue 
				if pvalue < 0.05:
					sig_pval.append(pos)


	return PTM_stats, sig_pval

def getBindingCore(options, refProt):

	# Import groundtruth and binding core 
	binders, bindingCores = importBindData(options)

	# Create array protein reference
	refProtStr = ''.join([refProt[AA] for AA in list(refProt.keys())])

	# Take binders with less than 10% of affinity, find the binding core,
	# and find indexes in protein of reference 
	strongBinders = list()
	coreIdxs = list()
	for binder in binders:
		if float(binders[binder]) <= 10 and binder in refProtStr:
			strongBinders.append(binder)
			core = bindingCores[binder]
			idx = re.search(core, refProtStr)
			coreIdxs.append(idx.span())
	
	return coreIdxs




def main():
    pass

if __name__ == "__main__":
    main()