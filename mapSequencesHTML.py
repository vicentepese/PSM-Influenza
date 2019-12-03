from Bio import Entrez
from Bio import SeqIO
import json
from subprocess import call
import csv
from collections import defaultdict, OrderedDict
import numpy as np 
import re
import os 
from operator import itemgetter
from markdown2 import Markdown
import random
random.seed(200)

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

def importData(options):

	# Import data  as list of lists
	data = list() 
	with open(options['files']['mergedFiltMASS'], 'r') as inFile:
		reader = csv.reader(inFile)
		for row in reader:
			data.append(row)
	
	return data

def parseSequences(options, data):

	# Initialize
	pos_range = options['pos_range']
	seqCount = defaultdict(lambda: defaultdict(int))
	seqInit = defaultdict(int)
	PTM_count = defaultdict(lambda: defaultdict(int))
	for seq in data:

		# Get initial position
		init_pos = int(seq[2])

		# Remove PTM and compute ending positon, else no PTM and compute ending position
		if '[' in seq[1]:
			
			# Get sequence and ending position
			AAseq = re.sub('\[.+?\]','',seq[1]).split('.')[1]
			end_pos = init_pos  + len(AAseq)-1
			AAseq = seq[1][2:-2]

			# Count number of PTMS and get index of PTM
			PTM_instances = re.findall('\[(.*?)\]', seq[1], re.DOTALL)
			for instance in PTM_instances:
				PTM_count[AAseq][instance] += 1
			
		else:
			AAseq = seq[1].split('.')[1]
			end_pos = init_pos + len(AAseq)-1

		# If initial position in the range, else check if any position overlap
		if not(end_pos < options['pos_range'][0]) and not(init_pos > options['pos_range'][1]):

			# Count and count in range
			seqCount[AAseq][seq[3]] += 1
			
			# Keep initial position 
			seqInit[AAseq] = init_pos

	return seqCount, seqInit, PTM_count

def findMutations(refProt, seq, init_pos):

	# Remove PTMs
	if '[' in seq:
		seq =  re.sub('\[.+?\]','',seq)

	# Create dictionnary 
	seqDict = {i: seq[i-init_pos] for i in range(init_pos,init_pos+len(seq))}

	# Check mutations 
	pos_mut = list()
	for pos in list(seqDict.keys()):
		if seqDict[pos] is not refProt[pos]:
			pos_mut.append(pos)
	
	# Convert to string index: zero-index base
	if len(pos_mut) > 0:
		pos_mut_idx = [pos-init_pos for pos in pos_mut]
	else:
		pos_mut_idx = []

	return pos_mut_idx

def getVaccineSample(seqCount):

	# Compute total count 
	vaccSample = defaultdict(int)
	for seq in list(seqCount.keys()):
		vaccSample['PAN'] += seqCount[seq]['PAN']
		vaccSample['ARP'] += seqCount[seq]['ARP']

	return vaccSample

def getRandomColor(options, PTM_count):

	# Get PTM types 
	uniquePTM = defaultdict(int)
	for seq in list(PTM_count.keys()):
		for ptm in list(PTM_count[seq].keys()):
			uniquePTM[ptm] += 1

	# Create color dictionnary	
	r = lambda: random.randint(75,200)
	color = {PTM: ['<span style=\"background: '+'#%02X%02X%02X; font-weight: bold' % (r(),r(),r()) + '\">',"</span>"] for PTM in list(uniquePTM.keys())}
	color['mutation'] = ['<span style=\"border-style: solid;  border-width: 2px\">','</span>']
	color['mutmatch'] = ['; border-style: solid; border-width: 2px']

	return color

def clearFragment(seqMark):
	if '[' in seqMark:
		seqMark = re.sub('\[.+?\]','',seqMark)
	return seqMark
			

def seqMutString(options, seq, pos_mut_idx, init_pos, seq_init_pos, seqCount, PTM_count, vaccSample, color):

	# Clear fragment 
	seqMark = clearFragment(seq)

	# Initialize 
	PTM = defaultdict(int)
	PTM_pos = list()
	PTM_type = list()

	# Get PTM type and location
	PTM_instances = re.findall('\[(.*?)\]', seq, re.DOTALL)
	PTM_idx = re.finditer('\[(.*?)\]', seq, re.DOTALL)
	idx_comp = 0
	for instance, idx in zip(PTM_instances, PTM_idx):
		PTM[instance] += 1
		PTM_pos.append(idx.start() -1 -idx_comp)
		PTM_type.append(instance)
		idx_comp += len(instance) +2

	# Remove negative and update number or ptms 
	PTM_pos = [item for item in PTM_pos if item >=0]

	# For each PTM, create Markdown string 
	PTM_pos_loop = PTM_pos
	for i in range(len(PTM_pos)):
		# Create string in Markdown and update indexing
		seqMark = seqMark[0:PTM_pos_loop[i]] + color[PTM_type[i]][0] + seqMark[PTM_pos_loop[i]] + color[PTM_type[i]][1] + seqMark[(PTM_pos_loop[i]+1):]
		PTM_pos_loop = [pos + len(color[PTM_type[i]][0]) + len(color[PTM_type[i]][1]) for pos in PTM_pos_loop]
	
	# Highlight mutations 
	if len(pos_mut_idx):
		if not any(mut in PTM_pos for mut in pos_mut_idx):
			pos_mut_idx  = [pos + (len(color[list(color.keys())[0]][0]) + len(color[list(color.keys())[0]][1]))*len(list(filter(lambda x: x < pos, PTM_pos))) for pos in pos_mut_idx]
			for i in range(0,len(pos_mut_idx)):
				seqMark = seqMark[0:pos_mut_idx[i]] + color['mutation'][0] + seqMark[pos_mut_idx[i]] + color['mutation'][1] + seqMark[(pos_mut_idx[i]+1):]
				pos_mut_idx = [pos + len(color['mutation'][0]) + len(color['mutation'][1]) for pos in pos_mut_idx]
		else:

			# Find matches of mutation and PTM (intersection) and remove them from original list
			pos_mut_match = list(set(pos_mut_idx).intersection(PTM_pos))
			for pos in pos_mut_match:
				pos_mut_idx.remove(pos)

			# Hightlight mutations
			pos_mut_idx  = [pos + (len(color[list(color.keys())[0]][0]) + len(color[list(color.keys())[0]][1]))*len(list(filter(lambda x: x < pos, PTM_pos))) for pos in pos_mut_idx]
			for i in range(0,len(pos_mut_idx)):
				seqMark = seqMark[0:pos_mut_idx[i]] + color['mutation'][0] + seqMark[pos_mut_idx[i]] + color['mutation'][1] + seqMark[(pos_mut_idx[i]+1):]
				pos_mut_idx = [pos + len(color['mutation'][0]) + len(color['mutation'][1]) for pos in pos_mut_idx]
			pos_mut_match  = [pos + (len(color[list(color.keys())[0]][0]) + len(color[list(color.keys())[0]][1])) * (len(list(filter(lambda x: x < pos, PTM_pos)))) + len(color[list(color.keys())[0]][0]) for pos in pos_mut_match]
			for i in range(0, len(pos_mut_match)):
				seqMark = seqMark[0:pos_mut_match[i]-2] + color['mutmatch'][0] + seqMark[pos_mut_match[i]-2:]
				pos_mut_match = [pos + len(color['mutmatch'][0]) for pos in pos_mut_match]
				
	# Append initial location and ARP and PAN proportion
	seqMark = '&nbsp;'*(np.absolute(init_pos-seq_init_pos)) + seqMark + ' \(ARP:{:.2%}' ' PAN:{:.2%}\)'.format(seqCount[seq]['ARP']/vaccSample['ARP'], seqCount[seq]['PAN']/vaccSample['PAN'])

	# Add PTMS at the end 
	for ptm in np.unique(np.asarray(PTM_type)):
		seqMark = seqMark + ' // ' + color[ptm][0] + ptm + color[ptm][1] 

	# Convert to HTML
	Markdowner = Markdown()
	seqHTML = Markdowner.convert(seqMark + '\n')

	return seqHTML


def mapOfSeqs(options, seqCount, seqInit, refProt, PTM_count, vaccSample, color):

	# Create array of positions: min initial pos to max ending pos
	init_pos = min(seqInit[seq] for seq in list(seqInit.keys()))

	# Indexed protein string 
	pos_range = options['pos_range']
	refProt_string = [refProt[pos] for pos in list(refProt.keys()) if (pos >= pos_range[0]) & (pos <= pos_range[1])]
	refProt_string = '&nbsp;'*(pos_range[0] - init_pos - 4) + str(pos_range[0]) +'.' + ''.join(refProt_string)  

	# Create string for each sequence 
	seqString = defaultdict(str)
	for seq, seq_init_pos in list(seqInit.items()):

		# Find mutations, get index for the string
		pos_mut_idx = findMutations(refProt, seq, seq_init_pos)

		# Create string in HTML format and store in dictionnary 
		seqHTML = seqMutString(options, seq, pos_mut_idx, init_pos, seq_init_pos, seqCount, PTM_count, vaccSample, color)
		seqString[seq] = seqHTML

	# Order by initial position
	seqInit = OrderedDict([(k, seqInit[k]) for k in sorted(seqInit, key=seqInit.get)])

	# Write 
	Markdowner = Markdown()
	with open(options['files']['seqMapOrex.html'],'w') as outFile:
		outFile.write('<style>' + '\n' + 'p { \n' + 'line-height:0.1; \n' + 'font-family: "Courier New", Courier, monospace; \n' + '} \n' + '\n' +"</style>" +'\n')
		outFile.write(Markdowner.convert(refProt_string + '\n'))
		outFile.write(Markdowner.convert('&nbsp; \n'))
		for seq in list(seqInit.keys()):
			outFile.write(seqString[seq])


def main():

	# Read options
	with open('options.json','r') as inFile:
		options = json.load(inFile)
	
	# Filter data 
	if len(os.listdir(options['folders']['FiltDataMASS'])):
		filterData(options)

	# Get reference protein 
	refProt = reference_retreive(options['refProt'])

	# Import MASSdata
	data = importData(options)

	# Map sequences and count
	seqCount, seqInit, PTM_count = parseSequences(options, data)

	# Write seqCount and seqCounRange
	with open(options['files']['seqCount'],'w') as outFile:
		writer = csv.writer(outFile)
		writer.writerow(['seq','count'])
		for seq in list(seqCount.keys()):
			writer.writerow([seq, seqCount[seq]])

	# Create random color 
	color = getRandomColor(options, PTM_count)

	# Compute vaccine sample 
	vaccSample = getVaccineSample(seqCount)

	# Create map of sequences 
	mapOfSeqs(options, seqCount, seqInit, refProt, PTM_count, vaccSample, color)



if __name__ == "__main__":
	main()
