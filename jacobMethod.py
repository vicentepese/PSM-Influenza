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
	handle = Entrez.efetch(db="protein", rettype="fasta",
						   retmode="text", id=proteinID)
	seq_record = SeqIO.read(handle, "fasta")

	# For each reference
	out_ref = OrderedDict()
	for i, j in enumerate(str(seq_record.seq)):
		out_ref[i+1] = j
	return out_ref


def importBindData(options):

	# Import binders
	binders = defaultdict(int)
	with open(options['files']['DQ0602Binders'], 'r') as inFile:
		next(inFile)
		reader = csv.reader(inFile)
		for row in reader:
			binders[row[2]] = float(row[3])

	# Import binding cores
	bindingCores = defaultdict(str)
	with open(options['files']['DQ0602BindingCores'], 'r') as inFile:
		next(inFile)
		reader = csv.reader(inFile, delimiter='\t')
		for row in reader:
			bindingCores[row[4]] = row[6]

	return binders, bindingCores


def getRandomColor(options, PTM_count):

	# Create color dictionnary
	r = lambda: random.randint(75, 200)
	color = {PTM: ['<span style=\"background: '+'#%02X%02X%02X; font-weight: bold' %
		(r(), r(), r()) + '\">', "</span>"] for PTM in list(PTM_count.keys())}
	color['ARP'] = ['<span style=\"color: #800000; font-weight: bold; \">', '</span>']
	color['PAN'] = [
		'<span style=\"color: #038f0a; font-weight: bold; \">', '</span>']
	color['bindingCore'] = [
		'<span style=\"background: green; font-weight: bold; \">', '</span>']

	return color


def div0(n, d):
	return n / d if d and n else 0


def countPositions(data):

	# Position count
	posCount = defaultdict(lambda: defaultdict(int))

	# Count number of occurences of an AA for each position
	for seq in data:

		# Clear data and count AA
		AAseq = re.sub('\[.+?\]', '', seq[1][2:-2])
		for i in range(0, len(AAseq)):
			posCount[int(seq[2])+i][AAseq[i]] += 1

	return posCount


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


def mapPTM(data, refProt, options):

	# Initialize
	seqPTM = defaultdict(lambda: defaultdict(
		lambda: defaultdict(int)))
	PTM_count = defaultdict(int)
	vaccSample = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

	for seq in data:

		# Initialize
		AAseq = seq[1][2:-2]
		AAnonPTM = re.sub('\[.+?\]', '', AAseq)
		init_pos = int(seq[2])

		# Count instances for each AA (position)
		for i in range(int(seq[2]), int(seq[2]) + len(AAnonPTM)):
			vaccSample[i][AAnonPTM[i-init_pos]][seq[3]] += 1

		# If there is a PTM
		if '[' in AAseq:

			# Get PTM indexes and types
			PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
			PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

			# For each PTM, find positions and append according to type
			idx_cumm = 0
			for instance, idx in zip(PTM_instances, PTM_idx):
				ptm_pos = init_pos + idx.start() - 1 - idx_cumm
				seqPTM[ptm_pos][instance][seq[3]] += 1
				PTM_count[instance] += 1
				idx_cumm += len(instance) + 2

	# Filter positions where there is no samples from any of the
	# vaccines
	seqPTM_loop = seqPTM
	for pos in list(seqPTM.keys()):
		for ptm in list(seqPTM[pos].keys()):
			if not(seqPTM[pos][ptm]['ARP']) or not(seqPTM[pos][ptm]['PAN']):
				del seqPTM[pos][ptm]
		if not(seqPTM[pos]):
			del seqPTM[pos]

	return seqPTM, vaccSample, PTM_count


def map2HTML(options, coreIdxs, refProt, PTM_stats, seqPTM, vaccSample, PTM_count):


	# Initialize
	PTM_HTML = list()
	markdowner = Markdown()
	color = getRandomColor(options, PTM_count)
	refProt = ''.join([refProt[pos] for pos in refProt])
	
	i = 0
	while i < len(refProt):

		# TODO: recheck this part, continue from here 
		# Create string of reference protein (taking 70 AA)
		refProtStr = refProt[i:i+70]
		count = 0 
		for core in coreIdxs:

			# If core is in that fragment of protein hightlight
			if core[0] in range(i, i + 70):
				# If no previous hightlight
				if count == 0:
					# Update core idxes
					core = [idx -i for idx in core]
					refProtStr = refProtStr[0:core[0]] + color['bindingCore'][0] + refProtStr[core[0]:core[1]] + \
						color['bindingCore'][1] + refProtStr[core[1]:]
					count += 1
				# If previous binding core in segment, update idx
				else:
					core = [idx - i +len(color['bindingCore'][0]) + len(color['bindingCore'][1]) for idx in core]
		
					refProtStr = refProtStr[0:core[0]] + color['bindingCore'][0] + refProtStr[core[0]:core[1]] + \
						color['bindingCore'][1] + refProtStr[core[1]:]
			elif core[1] in range(i, i + 70):
					# Update core idxes
					core = [idx -i for idx in core]
					core = [0 if idx < 0 else idx for idx in core]
					refProtStr = color['bindingCore'][0] + refProtStr[core[0]:core[1]] + \
						color['bindingCore'][1] + refProtStr[core[1]:]
					count += 1
			
		refProtStr = str(i+1) + '.' + '&nbsp;'*(6 -len(str(i))-1) + refProtStr + '\n'
		PTM_HTML.append(markdowner.convert(refProtStr))

		# Create ARP string
		# TODO: Continue from here
		ARP_str = color['ARP'][0] + 'ARP:&nbsp;&nbsp;' + color['ARP'][1]
		ARP_ptm = defaultdict(lambda: defaultdict(int))
		last_pos = 0
		for pos in range(i,i+70):
			if any(seqPTM[pos][ptm]['ARP'] for ptm in seqPTM[pos]):
				ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['ARP'][1] + refProt[pos-1]
				ARP_ptm[pos] = {ptm: seqPTM[pos][ptm] for ptm in seqPTM[pos] if seqPTM[pos][ptm]['ARP']}
				last_pos = pos - i
		ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(70 - last_pos) +  color['ARP'][1]
		PTM_HTML.append(markdowner.convert(ARP_str))

		# Create dictionnary of PTMS for that region of 70 
		ARP_ptm_str = '&nbsp;'*6
		while len(ARP_ptm) > 0:
			for pos in list(ARP_ptm.keys()):
				for ptm in list(ARP_ptm[pos].keys()):
					ARP_ptm_str  = ARP_ptm_str +  color[ptm][0] + ptm +  color['ARP'][1] + refProt[pos-1]



		# Create PAN string
		PAN_str = color['PAN'][0] + 'PAN:&nbsp;&nbsp;' + color['PAN'][1]
		last_pos = 0
		for pos in range(i,i+70):
			if any(seqPTM[pos][ptm]['PAN'] for ptm in seqPTM[pos]):
				PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['PAN'][1] + refProt[pos-1]
				last_pos = pos - i
		PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(70 - last_pos) +  color['PAN'][1]
		
		PTM_HTML.append(markdowner.convert(PAN_str))
		PTM_HTML.append(markdowner.convert('&nbsp;\n'))
		# Update index
		i += 70

	# Print and save
	with open(options['html']["scroll-template"], 'r') as inFile:
		with open(options['files']['mapJacob.html'], 'w') as outFile:
			for line in inFile:
				outFile.write(line)
			outFile.writelines(PTM_HTML)

def main():
	
	# Read options
	with open('options.json','r') as inFile:
		options = json.load(inFile)

	# Import data 
	data = importData(options)

	# Import protein of reference 
	refProt = reference_retreive(options['refProt'])

	# Get binding core and binding core positions
	coreIdxs = getBindingCore(options, refProt)

	# Get PTM positions, type and count 
	seqPTM, vaccSample, PTM_count = mapPTM(data, refProt, options)

	# Statistical test 
	PTM_stats = statisticalTest(options, seqPTM, vaccSample, refProt)

	# Create HTML output
	map2HTML(options, coreIdxs, refProt, PTM_stats, seqPTM, vaccSample, PTM_count)

	# Verbose
	print('stop')

if __name__ == "__main__":
	main()
