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


def getRandomColor(options):

	# Create color dictionnary
	color = {}
	color['ARP'] = ['<span style=\"color: #800000; font-weight: bold; \">', '</span>']
	color['PAN'] = ['<span style=\"color: #000782; font-weight: bold; \">', '</span>']
	color['strongBinder'] = ['<span style=\"background: #0F9D58; font-weight: bold; \">', '</span>']
	color['weakBinder'] = ['<span style=\"background: #F4B400; font-weight: bold; \">', '</span>']
	color['red'] = ['<span style=\"color: red; font-weight: bold; \">', '</span>']
	color['mut'] = ['<span style=\"font-weight: bold\">', '</span>']

	return color


def div0(n, d):
	return n / d if d and n else 0


def statisticalTest(options, seqMut, vaccSample, refProt):

	# Initialize 
	MUT_stats = defaultdict(lambda: defaultdict(lambda : defaultdict(int)))

	# For each position
	for pos in range(options['pos_range'][0], options['pos_range'][1]+1):

		if pos == 71:
			print('stop')
		# IF there is a mutation for both vaccines
		if len(list(seqMut[pos].keys())) >=1:
			for mut in list(seqMut[pos].keys()):
				if seqMut[pos][mut]['PAN'] and seqMut[pos][mut]['ARP']:

					# Create array 
					ptm_positive = [seqMut[pos][mut]['ARP'], seqMut[pos][mut]['PAN']]
					ptm_negative = [vaccSample[pos]['ARP'] - seqMut[pos][mut]['ARP'], \
						 vaccSample[pos]['PAN'] - seqMut[pos][mut]['PAN']]
					
					# Fisher test and append to output
					oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
					MUT_stats[pos][mut]['pvalue'] =  pvalue
					MUT_stats[pos][mut]['oddsratio'] = oddsratio

	return MUT_stats

def getBindingCore(options, refProt):

	# Import groundtruth and binding core
	binders, bindingCores = importBindData(options)

	# Create array protein reference
	refProtStr = ''.join([refProt[AA] for AA in list(refProt.keys())])

	# Take binders with less than 10% of affinity, find the binding core,
	# and find indexes in protein of reference
	coreIdxs = list()
	coreClass = list()
	for binder in binders:
		if float(binders[binder]) <= 20 and binder in refProtStr:
			core = bindingCores[binder]
			idx_binder = re.search(binder, refProtStr).span()
			idx_core = re.search(core, binder).span()
			idx =  [idx + idx_binder[0] for idx in idx_core]
			coreIdxs.append(idx)
			coreClass.append('strong')
		elif float(binders[binder]) <= 50 and binder in refProtStr:
			core = bindingCores[binder]
			idx_binder = re.search(binder, refProtStr).span()
			idx_core = re.search(core, binder).span()
			idx = [idx + idx_binder[0] for idx in idx_core]
			# Check for overlap
			if not any(idx[0] in range(coreRange[0], coreRange[1]) for coreRange in coreIdxs) and \
				not any(idx[1] in range(coreRange[0], coreRange[1]) for coreRange in coreIdxs):
				coreIdxs.append(idx)
				coreClass.append('weak')

	# Sort
	sortIdx = np.argsort([idx[0] for idx in coreIdxs])
	coreIdxs = [coreIdxs[idx] for idx in sortIdx]
	coreClass = [coreClass[idx] for idx in sortIdx]
	return coreIdxs, coreClass

def mapMutations(data, refProt, options):

	# Initialize
	seqMUT = defaultdict(lambda: defaultdict(lambda : defaultdict(int)))
	vaccSample = defaultdict(lambda: defaultdict((int)))

	for seq in data:

		# Initialize
		AAseq = seq[1][2:-2]
		AAnonPTM = re.sub('\[.+?\]', '', AAseq)
		init_pos = int(seq[2])

		# Check for mutationsseq[3]
		for AA, pos in zip(AAnonPTM, range(init_pos, init_pos + len(AAnonPTM))):

			# Count instances 
			vaccSample[pos][seq[3]] += 1

			# If there is a mutation append
			if AA is not refProt[pos]:
				seqMUT[pos][AA][seq[3]] += 1 
	
	# Filter positions where there is no samples from any of the
	# vaccines
	for pos in list(seqMUT.keys()):
		for mut in list(seqMUT[pos].keys()):
			if not(seqMUT[pos][mut]['ARP']) or not(seqMUT[pos][mut]['PAN']):
				del seqMUT[pos][mut]
		if len(seqMUT[pos]) < 1:
			del seqMUT[pos]

	return seqMUT, vaccSample


def map2HTML(options, coreIdxs, coreClass, refProt, MUT_stats, seqMut, vaccSample):


	# Initialize
	PTM_HTML = list()
	markdowner = Markdown()
	color = getRandomColor(options)
	refProt = ''.join([refProt[pos] for pos in refProt])
	
	i = 0
	while i < len(refProt):

		# Create string of reference protein (taking 70 AA)
		refProtStr = refProt[i:i+70]
		count = 0 
		for core, coreCl in zip(coreIdxs, coreClass):

			# If core is in that fragment of protein hightlight
			if core[0] in range(i, i + 70):
				# If no previous hightlight
				if count == 0:
					# Update core idxes
					core = [idx -i for idx in core]
					if coreCl == 'strong':
						refProtStr = refProtStr[0:core[0]] + color['strongBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['strongBinder'][1] + refProtStr[core[1]:]
						count += 1
					else:
						refProtStr = refProtStr[0:core[0]] + color['weakBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['weakBinder'][1] + refProtStr[core[1]:]
						count += 1
				# If previous binding core in segment, update idx
				else:
					if coreCl == 'strong':
						core = [idx - i + count*(len(color['strongBinder'][0]) + len(color['strongBinder'][1])) for idx in core]
						refProtStr = refProtStr[0:core[0]] + color['strongBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['strongBinder'][1] + refProtStr[core[1]:]
						count += 1
					else:
						core = [idx - i + count*(len(color['strongBinder'][0]) + len(color['strongBinder'][1])) for idx in core]
						refProtStr = refProtStr[0:core[0]] + color['weakBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['weakBinder'][1] + refProtStr[core[1]:]
						count += 1
			elif core[1] in range(i, i + 70):
					# Update core idxes
					core = [idx -i for idx in core]
					core = [0 if idx < 0 else idx for idx in core]
					if coreCl == 'strong':
						refProtStr = color['strongBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['strongBinder'][1] + refProtStr[core[1]:]
						count += 1
					else:
						refProtStr = color['weakBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['weakBinder'][1] + refProtStr[core[1]:]
						count += 1
		
		# Append
		refProtStr = str(i+1) + '.' + '&nbsp;'*(6 -len(str(i))-1) + refProtStr + '\n'
		PTM_HTML.append(markdowner.convert(refProtStr))

		# Create ARP string
		ARP_str = color['ARP'][0] + 'ARP:&nbsp;&nbsp;' + color['ARP'][1]
		ARP_mut = defaultdict(lambda: defaultdict(int))
		last_pos = 0
		for pos in range(i,i+70):
			if any(seqMut[pos][mut]['ARP'] for mut in list(seqMut[pos].keys())):
				ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['ARP'][1] + refProt[pos-1]
				ARP_mut[pos] = {mut: seqMut[pos][mut] for mut in seqMut[pos] if seqMut[pos][mut]['ARP']}
				last_pos = pos - i
		ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(70 - last_pos) +  color['ARP'][1]
		PTM_HTML.append(markdowner.convert(ARP_str))

		# Create PAN string
		PAN_str = color['PAN'][0] + 'PAN:&nbsp;&nbsp;' + color['PAN'][1]
		last_pos = 0
		for pos in range(i,i+70):
			if any(seqMut[pos][mut]['PAN'] for mut in seqMut[pos]):
				PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['PAN'][1] + refProt[pos-1]
				last_pos = pos - i
		PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(70 - last_pos) +  color['PAN'][1]
		
		PTM_HTML.append(markdowner.convert(PAN_str))

		# Create strings for each PTM positon and type 
		for pos in list(ARP_mut.keys()):
			for mut in list(ARP_mut[pos].keys()):
				ARP_prop = seqMut[pos][mut]['ARP']/vaccSample[pos]['ARP']
				ARP_samp = vaccSample[pos]['ARP']
				PAN_prop = seqMut[pos][mut]['ARP']/vaccSample[pos]['PAN']
				PAN_samp = vaccSample[pos]['PAN']
				ARP_ptm_str = '&nbsp;'*(pos -i -3+ 6) + \
						color['mut'][0] + mut +  color['mut'][1] + \
							'(ARP:{:.2%}({}), PAN:{:.2%}({}), '.format(ARP_prop, ARP_samp, PAN_prop, PAN_samp )
				if MUT_stats[pos][mut]['pvalue'] < 0.05:
					ARP_ptm_str = ARP_ptm_str + color['red'][0] + 'p={:.2}'.format(MUT_stats[pos][mut]['pvalue']) + '\n'
				else:
					ARP_ptm_str = ARP_ptm_str + 'p={:.2})'.format(MUT_stats[pos][mut]['pvalue']) + '\n'
				PTM_HTML.append(markdowner.convert(ARP_ptm_str))

		# Separate 
		PTM_HTML.append(markdowner.convert('&nbsp;\n'))

		# Update index
		i += 70

	# Print and save
	with open(options['html']["scroll-template"], 'r') as inFile:
		with open(options['files']['mutMapJacob.html'], 'w') as outFile:
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

	# Get binding cores and binding core positions
	coreIdxs, coreClass = getBindingCore(options, refProt)

	# Map mutations
	seqMut, vaccSample = mapMutations(data, refProt, options)

	# Compute Fisher exact test
	MUT_stats = statisticalTest(options, seqMut, vaccSample, refProt)

	# Create HTML output
	map2HTML(options, coreIdxs, coreClass, refProt, MUT_stats, seqMut, vaccSample)

if __name__ == "__main__":
	main()