import sys
import json
import re
import random
from markdown2 import Markdown
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from scipy import stats
from utils import getBindingCore, importBindData,\
	 importData, reference_retreive, div0, getBindingCore, getRandomColor

def statisticalTest(options, seqPTM, vaccSample, refProt):

	# Initialize 
	PTM_stats = defaultdict(lambda: defaultdict(lambda : defaultdict(lambda: defaultdict(int))))

	# For each position
	for pos in range(options['pos_range'][0], options['pos_range'][1]+1):

		# IF there is a ptm for both vaccines (ARP and FOC)
		if pos in list(seqPTM.keys()):
			for ptm in list(seqPTM[pos].keys()):
				if 'PAN' and 'ARP' in list(seqPTM[pos][ptm].keys()):

					# Create array 
					ptm_positive = [seqPTM[pos][ptm]['ARP'], seqPTM[pos][ptm]['PAN']]
					ptm_negative = [vaccSample[pos][refProt[pos]]['ARP'] - seqPTM[pos][ptm]['ARP'], \
						 vaccSample[pos][refProt[pos]]['PAN'] - seqPTM[pos][ptm]['PAN']]
					
					# Fisher test and append to output
					oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
					PTM_stats[pos][ptm]['ARP']['pvalue'] =  pvalue
					PTM_stats[pos][ptm]['ARP']['oddsratio'] = oddsratio

				if 'PAN' and 'FOC' in list(seqPTM[pos][ptm].keys()):

					# Create array 
					ptm_positive = [seqPTM[pos][ptm]['FOC'], seqPTM[pos][ptm]['PAN']]
					ptm_negative = [vaccSample[pos][refProt[pos]]['FOC'] - seqPTM[pos][ptm]['FOC'], \
						 vaccSample[pos][refProt[pos]]['PAN'] - seqPTM[pos][ptm]['PAN']]
					
					# Fisher test and append to output
					oddsratio, pvalue = stats.fisher_exact([ptm_positive, ptm_negative])
					PTM_stats[pos][ptm]['FOC']['pvalue'] =  pvalue
					PTM_stats[pos][ptm]['FOC']['oddsratio'] = oddsratio

	return PTM_stats

def mapPTM(data, refProt, options):

	# Initialize
	seqPTM = defaultdict(lambda: defaultdict(
		lambda: defaultdict(int)))
	PTM_count = defaultdict(lambda: defaultdict(int))
	vaccSample = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

	# For each sequence
	for seq in data:

		# Initialize: sequence with and without PTM, initial position
		AAseq = seq[1][2:-2]
		AAnonPTM = re.sub('\[.+?\]', '', AAseq)
		init_pos = int(seq[2])

		# Count instances for each AA (position)
		for i in range(init_pos, init_pos + len(AAnonPTM)):
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
				PTM_count[instance][seq[3]] += 1
				idx_cumm += len(instance) + 2

	# Filter positions where there is no samples from any of the vaccines
	for pos in list(seqPTM.keys()):
		for ptm in list(seqPTM[pos].keys()):
			if not(seqPTM[pos][ptm]['ARP'] and seqPTM[pos][ptm]['PAN']) \
				and not(seqPTM[pos][ptm]['FOC'] and seqPTM[pos][ptm]['PAN']):
				del seqPTM[pos][ptm]
		if len(seqPTM[pos]) < 1:
			del seqPTM[pos]

	return seqPTM, vaccSample, PTM_count


def map2HTML(options, coreIdxs, coreClass, refProt, PTM_stats, seqPTM, vaccSample, PTM_count):


	# Initialize
	PTM_HTML = list()
	markdowner = Markdown()
	color = getRandomColor(options, PTM_count=PTM_count)
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

				# If no previous core
				if count == 0:

					# Update core idxes, and highlight based on class
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
			if pos in list(seqPTM.keys()):
				if any(seqPTM[pos][ptm]['PAN'] for ptm in seqPTM[pos]):
					PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['PAN'][1] + refProt[pos-1]
					last_pos = pos - i
		PAN_str  = PAN_str +  color['PAN'][0] + '&mdash;'*(70 - last_pos) +  color['PAN'][1]
		PTM_HTML.append(markdowner.convert(PAN_str))

		# Create ARP string, highlighting positions of PTMs, and append
		ARP_str = color['ARP'][0] + 'ARP:&nbsp;&nbsp;' + color['ARP'][1]
		ptm_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
		last_pos = 0
		for pos in range(i,i+70):
			if pos in list(seqPTM.keys()):
				if any(seqPTM[pos][ptm]['ARP'] for ptm in seqPTM[pos]):
					ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['ARP'][1] + refProt[pos-1]
					for ptm in seqPTM[pos]:
						ptm_dict[pos][ptm]['ARP'] = seqPTM[pos][ptm]['ARP']
						last_pos = pos - i
		ARP_str  = ARP_str +  color['ARP'][0] + '&mdash;'*(70 - last_pos) +  color['ARP'][1]
		PTM_HTML.append(markdowner.convert(ARP_str))

		# Create FOC string, highlighting positions of PTMs, and append
		FOC_str = color['FOC'][0] + 'FOC:&nbsp;&nbsp;' + color['FOC'][1]
		last_pos = 0
		for pos in range(i,i+70):
			if pos in list(seqPTM.keys()):
				if any(seqPTM[pos][ptm]['FOC'] for ptm in seqPTM[pos]):
					FOC_str  = FOC_str +  color['FOC'][0] + '&mdash;'*(pos - last_pos -1 - i) +  color['FOC'][1] + refProt[pos-1]
					for ptm in seqPTM[pos]:
						ptm_dict[pos][ptm]['FOC'] = seqPTM[pos][ptm]['FOC'] 
						last_pos = pos - i
		FOC_str  = FOC_str +  color['FOC'][0] + '&mdash;'*(70 - last_pos) +  color['FOC'][1]
		PTM_HTML.append(markdowner.convert(FOC_str))

		# Create strings for each PTM positon and type 
		for pos in list(ptm_dict.keys()):
			for ptm in list(ptm_dict[pos].keys()):
				for vacc in list(ptm_dict[pos][ptm].keys()):
					vacc_prop = seqPTM[pos][ptm][vacc]/vaccSample[pos][refProt[pos-1]][vacc]
					vacc_samp = vaccSample[pos][refProt[pos-1]][vacc]
					PAN_prop = seqPTM[pos][ptm]['PAN']/vaccSample[pos][refProt[pos-1]]['PAN'] 
					PAN_samp = vaccSample[pos][refProt[pos-1]]['PAN']
					PAN_ptm_str = '&nbsp;'*(pos -i -3+ 6) + \
							color[ptm][0] + ptm +  color[ptm][1] + \
								'(' + vacc + ':{:.2%}({}),PAN:{:.2%}({}),'.format(vacc_prop, vacc_samp, PAN_prop, PAN_samp)
					if pos in list(PTM_stats.keys()) and vacc in list(PTM_stats[pos][ptm].keys()) \
						and PTM_stats[pos][ptm][vacc]['pvalue'] < 0.05:
						PAN_ptm_str = PAN_ptm_str + color['red'][0] + 'p={:.2}'.format(PTM_stats[pos][ptm][vacc]['pvalue']) + '\n'
					elif pos in list(PTM_stats.keys()) and vacc in list(PTM_stats[pos][ptm].keys()):
						PAN_ptm_str = PAN_ptm_str + 'p={:.2})'.format(PTM_stats[pos][ptm][vacc]['pvalue']) + '\n'
					PTM_HTML.append(markdowner.convert(PAN_ptm_str))

		# Separate 
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

	# Jacob method is designed to print the entire sequence of the protein reference 
	options['pos_range'] = [1, 566]

	# Import data 
	data = importData(options)

	# Import protein of reference 
	refProt = reference_retreive(options['refProt'])

	# Get binding core and binding core positions
	coreIdxs, coreClass = getBindingCore(options, refProt)

	# Get PTM positions, type and count 
	seqPTM, vaccSample, PTM_count = mapPTM(data, refProt, options)

	# Statistical test 
	PTM_stats = statisticalTest(options, seqPTM, vaccSample, refProt)

	# Create HTML output
	map2HTML(options, coreIdxs, coreClass, refProt, PTM_stats, seqPTM, vaccSample, PTM_count)

if __name__ == "__main__":
	main()
