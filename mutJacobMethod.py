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
from utils import getBindingCore, importBindData,\
	 importData, reference_retreive, div0, getBindingCore, getRandomColor

def statisticalTest(options, seqMut, vaccSample, refProt):

	# Initialize 
	MUT_stats = defaultdict(lambda: defaultdict(lambda : defaultdict(lambda: defaultdict(int))))

	# For each position
	for pos in range(options['pos_range'][0], options['pos_range'][1]+1):

		if pos in list(seqMut.keys()):
			for ptm in list(seqMut[pos].keys()):
				if 'PAN' and 'ARP' in list(seqMut[pos][ptm].keys()):

					# Create array 
					ptm_positive = [seqMut[pos][ptm]['ARP'], seqMut[pos][ptm]['PAN']]
					ptm_negative = [vaccSample[pos]['ARP'] - seqMut[pos][ptm]['ARP'], \
						 vaccSample[pos]['PAN'] - seqMut[pos][ptm]['PAN']]
					
					# Fisher test and append to output
					oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
					MUT_stats[pos][ptm]['ARP']['pvalue'] =  pvalue
					MUT_stats[pos][ptm]['ARP']['oddsratio'] = oddsratio

				if 'PAN' and 'FOC' in list(seqMut[pos][ptm].keys()):

					# Create array 
					ptm_positive = [seqMut[pos][ptm]['FOC'], seqMut[pos][ptm]['PAN']]
					ptm_negative = [vaccSample[pos]['FOC'] - seqMut[pos][ptm]['FOC'], \
						 vaccSample[pos]['PAN'] - seqMut[pos][ptm]['PAN']]
					
					# Fisher test and append to output
					oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
					MUT_stats[pos][ptm]['FOC']['pvalue'] =  pvalue
					MUT_stats[pos][ptm]['FOC']['oddsratio'] = oddsratio

	return MUT_stats


def mapMutations(data, refProt, options):

	# Initialize outputs
	seqMUT = defaultdict(lambda: defaultdict(lambda : defaultdict(int)))
	vaccSample = defaultdict(lambda: defaultdict((int)))

	# For each sequence
	for seq in data:

		# Initialize: sequence with and without PTM, initial position
		AAseq = seq[1][2:-2]
		AAnonPTM = re.sub('\[.+?\]', '', AAseq)
		init_pos = int(seq[2])

		# Check for mutations
		for AA, pos in zip(AAnonPTM, range(init_pos, init_pos + len(AAnonPTM))):

			# Count instances 
			vaccSample[pos][seq[3]] += 1

			# If there is a mutation append
			if AA is not refProt[pos]:
				seqMUT[pos][AA][seq[3]] += 1 
	
	# Filter positions where there is no samples from any of the
	# vaccines
	for pos in list(seqMUT.keys()):
		for ptm in list(seqMUT[pos].keys()):
			if not(seqMUT[pos][ptm]['ARP'] and seqMUT[pos][ptm]['PAN']) \
				and not(seqMUT[pos][ptm]['FOC'] and seqMUT[pos][ptm]['PAN']):
				del seqMUT[pos][ptm]
		if len(seqMUT[pos]) < 1:
			del seqMUT[pos]

	return seqMUT, vaccSample


def map2HTML(options, coreIdxs, coreClass, refProt, MUT_stats, seqMut, vaccSample):

	# Initialize
	PTM_HTML = list()
	markdowner = Markdown()
	color = getRandomColor(options)
	refProt = ''.join([refProt[pos] for pos in refProt])
	
	# In blocks of 70, while smaller than the length of the protein of reference
	i = 0
	while i < len(refProt):

		# Create string of reference protein (taking 70 AA)
		refProtStr = refProt[i:i+70]
		count = 0 

		# For each binding core and class
		for core, coreCl in zip(coreIdxs, coreClass):

			# If initial position of the core overlaps with that fragment
			if core[0] in range(i, i + 70):

				# If no previous hightlight
				if count == 0:

					#  Update core idxes, and highlight based on class
					core = [idx -i for idx in core]
					if coreCl == 'strong':
						refProtStr = refProtStr[0:core[0]] + color['strongBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['strongBinder'][1] + refProtStr[core[1]:]
						count += 1
					else:
						refProtStr = refProtStr[0:core[0]] + color['weakBinder'][0] + refProtStr[core[0]:core[1]] + \
							color['weakBinder'][1] + refProtStr[core[1]:]
						count += 1

				# If previous binding core in segment, update idx and highlight based on class
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

			# If ending position of the core overlaps with the fragment: same	
			elif core[1] in range(i, i + 70):
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
		
		# Append to HTML output
		refProtStr = str(i+1) + '.' + '&nbsp;'*(6 -len(str(i))-1) + refProtStr + '\n'
		PTM_HTML.append(markdowner.convert(refProtStr))

		# Create PAN string: same as ARP string
		PAN_str = color['PAN'][0] + 'PAN:&nbsp;&nbsp;' + color['PAN'][1]
		last_pos = 0
		for pos in range(i,i+70):
			if pos in list(seqMut.keys()):
				if any(seqMut[pos][mut]['PAN'] for mut in seqMut[pos]):
					PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['PAN'][1] + refProt[pos-1]
					last_pos = pos - i
		PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(70 - last_pos) +  color['PAN'][1]
		PTM_HTML.append(markdowner.convert(PAN_str))

		# Create ARP string, highlighting positions of PTMs, and append
		ARP_str = color['ARP'][0] + 'ARP:&nbsp;&nbsp;' + color['ARP'][1]
		mut_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
		last_pos = 0
		for pos in range(i,i+70):
			if pos in list(seqMut.keys()):
				if any(seqMut[pos][mut]['ARP'] for mut in seqMut[pos]):
					ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['ARP'][1] + refProt[pos-1]
					for mut in seqMut[pos]:
						mut_dict[pos][mut]['ARP'] = seqMut[pos][mut]['ARP']
						last_pos = pos - i
		ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(70 - last_pos) +  color['ARP'][1]
		PTM_HTML.append(markdowner.convert(ARP_str))

		# Create FOC string, highlighting positions of PTMs, and append
		FOC_str = color['FOC'][0] + 'FOC:&nbsp;&nbsp;' + color['FOC'][1]
		last_pos = 0
		for pos in range(i,i+70):
			if pos in list(seqMut.keys()):
				if any(seqMut[pos][mut]['FOC'] for mut in seqMut[pos]):
					FOC_str  = FOC_str +  color['FOC'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['FOC'][1] + refProt[pos-1]
					for mut in seqMut[pos]:
						mut_dict[pos][mut]['FOC'] = seqMut[pos][mut]['FOC']
						last_pos = pos - i
		FOC_str  = FOC_str +  color['FOC'][0] + '&mdash;'*(70 - last_pos) +  color['FOC'][1]
		PTM_HTML.append(markdowner.convert(FOC_str))

		# Create strings for each PTM positon and type 
		for pos in list(mut_dict.keys()):
			for mut in list(mut_dict[pos].keys()):
				for vacc in list(mut_dict[pos][mut].keys()):
					if mut_dict[pos][mut][vacc] > 0:
						vacc_prop = seqMut[pos][mut][vacc]/vaccSample[pos][vacc]
						vacc_samp = vaccSample[pos][vacc]
						PAN_prop = seqMut[pos][mut]['PAN']/vaccSample[pos]['PAN'] 
						PAN_samp = vaccSample[pos]['PAN']
						PAN_mut_str = '&nbsp;'*(pos -i -3+ 6) + \
								color['mut'][0] + mut +  color['mut'][1] + \
									'(' + vacc + ':{:.2%}({}),PAN:{:.2%}({}),'.format(vacc_prop, vacc_samp, PAN_prop, PAN_samp)
						if pos in list(MUT_stats.keys()) and vacc in list(MUT_stats[pos][mut].keys()) \
							and MUT_stats[pos][mut][vacc]['pvalue'] < 0.05:
							PAN_mut_str = PAN_mut_str + color['red'][0] + 'p={:.2}'.format(MUT_stats[pos][mut][vacc]['pvalue']) + '\n'
						elif pos in list(MUT_stats.keys()) and vacc in list(MUT_stats[pos][mut].keys()):
							PAN_mut_str = PAN_mut_str + 'p={:.2})'.format(MUT_stats[pos][mut][vacc]['pvalue']) + '\n'
						PTM_HTML.append(markdowner.convert(PAN_mut_str))

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