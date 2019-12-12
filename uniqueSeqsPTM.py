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

def checkOverlap(seq_init_pos, seq_end_pos, coreIdxs):
	# If there is a partial overlap of at least one AA
	for core in coreIdxs:
		if seq_init_pos in range(core[0], core[1]+1) or seq_end_pos in range(core[0], core[1]+1):
			return True, core

	return False, None

def mapSeqPTM(data, refprot, options):

	# Initialize
	seqPTM = defaultdict(lambda: defaultdict(lambda : defaultdict(lambda : defaultdict(int))))
	seqCount = defaultdict(lambda: defaultdict(int))
	seqInit = defaultdict(int)
	PTM_count = defaultdict(int)

	# For each sequence get sequence of AA (without PTM), initial and ending position
	for seq in data:

		AAseq = re.sub('\[.+?\]','',seq[1])[2:-2]
		init_pos = int(seq[2])
		end_pos = init_pos + len(AAseq)

		# If sequence overlaps with range
		if not(end_pos < options['pos_range'][0]) and not(init_pos > options['pos_range'][1]):
			seqCount[AAseq][seq[3]] += 1
			seqInit[AAseq] = int(seq[2])

			# If sequence has PTM, save instances and indexes, else save sequence
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

def seq2HTML(options, seqPTM, seqCount, seqInit, PTM_count, refProt, coreIdxs):

	# Initialize 
	init_pos = min(seqInit[seq] for seq in list(seqInit.keys()))
	seqHTML = defaultdict(str)
	Markdowner = Markdown()

	# Get random colors
	color = getRandomColor(options, PTM_count)

	# For each sequence
	for seq in list(seqPTM.keys()):

		# Initialize
		seqMark = seq
		PTM_pos_loop = list(seqPTM[seq].keys())
		PTM_pos_loop.sort()
		seqPos = list(seqPTM[seq].keys())
		seqPos.sort()
		seq_init_pos = seqInit[seq]
		seq_end_pos = seq_init_pos + len(seq)

		# Compute Fisher extact test for each PTM between vaccines
		PTM_stats, sig_pval = statisticalTest(options, seqPTM[seq], seqCount, seq)

		# Highlight binding cores:
		check, core = checkOverlap(seq_init_pos, seq_end_pos, coreIdxs)
		if check:
			core = [idx + 1 - seq_init_pos for idx in core]
			core = [0 if pos < 0 else pos for pos in core]
			seqMark = seqMark[0:core[0]] + color['bindingCore'][0] + seqMark[core[0]:core[1]] + \
					color['bindingCore'][1] + seqMark[core[1]:]
			for i in range(len(PTM_pos_loop)):
				if PTM_pos_loop[i] in range(core[0], core[1]):
					PTM_pos_loop[i] = PTM_pos_loop[i] + len(color['bindingCore'][0])
				elif PTM_pos_loop[i] < core[0]:
					continue
				else:
					PTM_pos_loop[i] =PTM_pos_loop[i] = PTM_pos_loop[i] + len(color['bindingCore'][0]) + len(color['bindingCore'][1])

		# Create markdown string (highlight PTMs in the sequence)
		for i in range(0,len(PTM_pos_loop)):
			if seqPos[i] in sig_pval:
				seqMark = seqMark[0:PTM_pos_loop[i]] + color['red'][0] + seqMark[PTM_pos_loop[i]] + \
					color['red'][1] + seqMark[(PTM_pos_loop[i]+1):]
				PTM_pos_loop = [pos + len(color['red'][0]) + len(color['red'][1]) for pos in PTM_pos_loop]
			else:
				seqMark = seqMark[0:PTM_pos_loop[i]] + color['orange'][0] + seqMark[PTM_pos_loop[i]] + \
					color['orange'][1] + seqMark[(PTM_pos_loop[i]+1):]
				PTM_pos_loop = [pos + len(color['orange'][0]) + len(color['orange'][1]) for pos in PTM_pos_loop]

		# TODO: optimize this if-else statement
		# Append initial location and ARP and PAN proportion 
		if options['pos_range'][0] > 1:
			seqMark = '&nbsp;'*(np.absolute(init_pos-seq_init_pos)) + seqMark + \
				'(ARP: {0}, PAN: {1}'.format(seqCount[seq]['ARP'], seqCount[seq]['PAN']) + ')'
		else:
			seqMark = '&nbsp;'*(np.absolute(init_pos-seq_init_pos) - len(str(seq_init_pos)) - 1) + str(seq_init_pos) + '.' + seqMark + \
				'(ARP: {0}, PAN: {1}'.format(seqCount[seq]['ARP'], seqCount[seq]['PAN']) + ')'

		for pos in seqPos:
			seqMark = seqMark +  ' // ' + '__' + str(pos) + '__: '
			for ptm in list(seqPTM[seq][pos].keys()):
				if PTM_stats[pos][ptm]:
					seqMark = seqMark + color[ptm][0] + ptm + color[ptm][1] + \
						' \(ARP:{:.2%}' ' PAN:{:.2%}, p = {:.2}\) '.format(div0(seqPTM[seq][pos][ptm]['ARP'],seqCount[seq]['ARP']),\
							div0(seqPTM[seq][pos][ptm]['PAN'],seqCount[seq]['PAN']), PTM_stats[pos][ptm]['pvalue'])
				else:
					seqMark = seqMark + color[ptm][0] + ptm + color[ptm][1] + \
						' \(ARP:{:.2%}' ' PAN:{:.2%}\) '.format(div0(seqPTM[seq][pos][ptm]['ARP'],seqCount[seq]['ARP']),\
							div0(seqPTM[seq][pos][ptm]['PAN'],seqCount[seq]['PAN']))

		seqHTML[seq] = Markdowner.convert(seqMark + '\n')
	
	# Sort initial position
	seqInit = OrderedDict([(k, seqInit[k]) for k in sorted(seqInit, key=seqInit.get)])

	# String for protein of reference 
	pos_range = options['pos_range']
	refProt_string = [refProt[pos] for pos in list(refProt.keys()) if (pos >= pos_range[0]) & (pos <= pos_range[1])]
	refProt_string = '&nbsp;'*(pos_range[0] - init_pos - 4) + str(pos_range[0]) +'.' + ''.join(refProt_string)

	# Save 	
	with open(options['html']["scroll-template"], 'r') as inFile:	
		with open(options['files']['uniqueSeqsPTM.html'],'w') as outFile:
			for line in inFile:
				outFile.write(line)
			refProt_string = Markdowner.convert(refProt_string + '\n')
			refProt_string = refProt_string
			outFile.write(Markdowner.convert(refProt_string + '\n'))
			outFile.write(Markdowner.convert('&nbsp; \n'))
			for seq in seqInit.keys():
				if seqHTML[seq]:
					seq_string = seqHTML[seq]
					outFile.write(seq_string)
				else:
					continue
			

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

	# Get binding core and binding core positions
	coreIdxs = getBindingCore(options, refProt)

	# Compute HTML document 
	seq2HTML(options, seqPTM, seqCount, seqInit, PTM_count, refProt, coreIdxs)

if __name__ == "__main__":
	main()