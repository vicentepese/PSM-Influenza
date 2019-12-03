import numpy as np 
import sys 
import os 
import json 
import csv 
import re
import random 
from markdown2 import Markdown
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict, OrderedDict

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

def getRandomColor(options, PTM_count):

	# Get PTM types 


	# Create color dictionnary	
	r = lambda: random.randint(75,200)
	color = {PTM: ['<span style=\"background: '+'#%02X%02X%02X; font-weight: bold' % (r(),r(),r()) + '\">',"</span>"] for PTM in list(PTM_count.keys())}
	color['mutation'] = ['<span style=\"border-style: solid;  border-width: 2px\">','</span>']
	color['mutmatch'] = ['; border-style: solid; border-width: 2px']
	color['red'] = ['<span style=\"background: red; font-weight: bold\"', '</span>']

	return color

def mapSeqPTM(data, refprot, options):

	# Initialize
	seqPTM = defaultdict(lambda: defaultdict(lambda : defaultdict(lambda : defaultdict(int))))
	seqCount = defaultdict(lambda: defaultdict(int))
	seqInit = defaultdict(int)
	PTM_count = defaultdict(int)

	# For each sequence get sequence of AA (without PTM), initial and ending position
	for seq in data:
		AAseq = re.sub('\[.+?\]','',seq[1])[2:-2]
		if 'AFAMERNAGSGF' in AAseq:
			print('stop')
		init_pos = int(seq[2])
		end_pos = init_pos + len(AAseq)

		# If sequence overlaps with range
		if not(end_pos < options['pos_range'][0]) and not(init_pos > options['pos_range'][1]):
			seqCount[AAseq][seq[3]] += 1
			seqInit[AAseq] = int(seq[2])

			# If sequence has PTM, save instances and indexes
			if '[' in seq[1][2:-2]:
				PTM_instances = re.findall('\[(.*?)\]', seq[1][2:-2], re.DOTALL)
				PTM_idx = re.finditer('\[(.*?)\]', seq[1][2:-2], re.DOTALL)
				idx_cumm = 0

				for instance, idx in zip(PTM_instances, PTM_idx):
					if idx.start() > 0:
						seqPTM[AAseq][idx.start() - 1 - idx_cumm][instance][seq[3]] += 1
						PTM_count[instance] += 1
					idx_cumm += len(instance) + 2

	return seqPTM, seqCount, seqInit, PTM_count

def seq2HTML(options, seqPTM, seqCount, seqInit, PTM_count, refProt):

	# Initialize 
	init_pos = min(seqInit[seq] for seq in list(seqInit.keys()))
	seqHTML = list()
	Markdowner = Markdown()
	# Get random colors
	color = getRandomColor(options, PTM_count)

	# For each sequence
	for seq in list(seqPTM.keys()):
		seqMark = seq
		PTM_pos_loop = list(seqPTM[seq].keys())
		for pos in PTM_pos_loop:
			for ptm in list(seqPTM[seq][pos].keys()):
				seqMark = seqMark[0:pos] + color['red'][0] + seqMark[pos] + color['red'][1] + seqMark[(pos+1):]
			PTM_pos_loop = [pos + len(color[ptm][0]) + len(color[ptm][1]) for pos in PTM_pos_loop]

		# Append initial location and ARP and PAN proportion 
		seq_init_pos = seqInit[seq]
		seqMark = '&nbsp;'*(np.absolute(init_pos-seq_init_pos)) + seqMark 
		for pos in list(seqPTM[seq].keys()):
			for ptm in list(seqPTM[seq][pos].keys()):
				seqMark = seqMark + '//' + '__' + str(pos) + '__: ' + color[ptm][0] + ptm + color[ptm][1] + ' \(ARP:{:.2%}' ' PAN:{:.2%}\)'.format(seqPTM[seq][pos][ptm]['ARP']/PTM_count['ARP'], seqPTM[seq][pos][ptm]['PAN']/PTM_count['PAN'])
		seqHTML.append(Markdowner.convert(seqMark + '\n'))

	# String for protein of reference 
	pos_range = options['pos_range']
	refProt_string = [refProt[pos] for pos in list(refProt.keys()) if (pos >= pos_range[0]) & (pos <= pos_range[1])]
	refProt_string = '&nbsp;'*(pos_range[0] - init_pos - 4) + str(pos_range[0]) +'.' + ''.join(refProt_string)


	# Save 	
	with open('test.html','w') as outFile:
		outFile.write('<style>' + '\n' + 'p { \n' + 'line-height:0.1; \n' + 'font-family: "Courier New", Courier, monospace; \n' + '} \n' + '\n' +"</style>" +'\n')
		outFile.write(Markdowner.convert(refProt_string + '\n'))
		outFile.write(Markdowner.convert('&nbsp; \n'))
		outFile.writelines(seqHTML)

def main():

	# Read options 
	with open('options.json', 'r') as inFile:
		options = json.load(inFile)

	# Import  data
	data = importData(options)

	# Import protein of reference 
	refProt = reference_retreive(options['refProt'])

	# Count PTMs for each unique sequence
	seqPTM, seqCount, seqInit, PTM_count = mapSeqPTM(data, refProt, options)

	# Compute HTML document 
	seq2HTML(options, seqPTM, seqCount, seqInit, PTM_count, refProt)

if __name__ == "__main__":
	main()