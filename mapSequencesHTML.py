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
	with open(options['files']['mergedFiltMASS'], 'r') as inFile:
		reader = csv.reader(inFile)
		for row in reader:
			AFLSA.append(row)
	
	return AFLSA

def mapSeqs(options, AFLSA):

	# For each sequence 
	pos_range = options['pos_range']

	seqMap = defaultdict(list)
	for seq in AFLSA:

		# Get initial position
		init_pos = int(seq[2])

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

				seqMap[str(pos)].append(AA)
		
		elif end_pos >= pos_range[0] and end_pos <= pos_range[1]:

			# Map each position 
			for pos, AA in zip(range(init_pos, (end_pos + 1)), AAseq):

				seqMap[str(pos)].append(AA)
			
	return seqMap

def countSequences(options, AFLSA):

	# For each sequence 
	pos_range = options['pos_range']
	seqCount = defaultdict(lambda: defaultdict(int))
	seqCountRange = defaultdict(lambda: defaultdict(int))
	seqInit = defaultdict(int)
	for seq in AFLSA:

		 # Get initial position
		init_pos = int(seq[2])

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

			# Count 
			seqCount[AAseq][seq[3]] += 1

			# Count in range 
			seqCountRange[AAseq[0:(pos_range[1]-init_pos)]][seq[3]] += 1

			# Keep initial position 
			seqInit[AAseq] = init_pos
		
		elif end_pos >= pos_range[0] and end_pos <= pos_range[1]:

			# Count 
			seqCount[AAseq][seq[3]] += 1

			# Count in range 
			seqCountRange[AAseq[pos_range[0] - init_pos:]][seq[3]] += 1

			# Keep initial position 
			seqInit[AAseq] = init_pos

	return seqCount, seqCountRange, seqInit

def findMutations(refProt, seq, init_pos):

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

def seqMutString(options, seq, pos_mut_idx, init_pos, seq_init_pos, seqCount):

	# For each mutation, add two _ before and after the idx
	seqMark = seq 
	for i in range(0, len(pos_mut_idx)):

		# Create string in Markdown 
		seqMark = seqMark[0:(pos_mut_idx[i])] + options['MarkdownHigh']['init'] + seqMark[pos_mut_idx[i]] + options['MarkdownHigh']['end'] + seqMark[(pos_mut_idx[i]+1):]
		
		# Update indexing
		pos_mut_idx = [pos + len(options['MarkdownHigh']['init']) + len(options['MarkdownHigh']['end']) for pos in pos_mut_idx]

	# Add spaces until initial position 
	seqMark = '&nbsp;'*(np.absolute(init_pos-seq_init_pos)) + seqMark + '(ARP:' + str(seqCount[seq]['ARP']) + ', PAN:' + str(seqCount[seq]['PAN']) +  ')'
	# Converto HMLT 
	Markdowner = Markdown()
	seqHTML = Markdowner.convert(seqMark + '\n')

	return seqHTML


def mapOfSeqs(options, seqCount, seqInit, refProt):

	# Create array of positions: min initial pos to max ending pos
	init_pos = min(seqInit[seq] for seq in list(seqInit.keys()))
	ending_pos = max(seqInit[seq] + len(seq) for seq in list(seqInit.keys()))

	# Indexed protein string 
	pos_range = options['pos_range']
	refProt_string = [refProt[pos] for pos in list(refProt.keys()) if (pos >= pos_range[0]) & (pos <= pos_range[1])]
	refProt_string = '&nbsp;'*(pos_range[0] - init_pos) + ''.join(refProt_string) + ' '*(ending_pos - pos_range[1])

	# Create string for each sequence 
	seqString = defaultdict(str)
	for seq, seq_init_pos in list(seqInit.items()):

		# Find mutations, get index for the string
		pos_mut_idx = findMutations(refProt, seq, seq_init_pos)

		# Create string in HTML format and store in dictionnary 
		seqHTML = seqMutString(options, seq, pos_mut_idx, init_pos, seq_init_pos, seqCount)
		seqString[seq] = seqHTML


	# Order by initial position
	seqInit = OrderedDict([(k, seqInit[k]) for k in sorted(seqInit, key=seqInit.get)])

	# Write 
	Markdowner = Markdown()
	with open(options['files']['seqMap.html'],'w') as outFile:
		outFile.write('<style>' + '\n' + 'p { \n' + 'line-height:0.1; \n' + 'font-family: "Courier New", Courier, monospace; \n' + '} \n' + '\n' +"</style>" +'\n')
		outFile.write(Markdowner.convert(refProt_string + '\n'))
		outFile.write(Markdowner.convert('\n'))
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
	AFLSA = importAFLSA(options)

	# Map sequences 
	seqMap = mapSeqs(options, AFLSA)

	for pos in list(seqMap.keys()):
		seqMap[pos] = list(np.unique(np.asarray(seqMap[pos])))

	# Count sequences
	seqCount, seqCountRange, seqInit = countSequences(options, AFLSA)

	# Write seqCount and seqCounRange
	with open(options['files']['seqCount'],'w') as outFile:
		writer = csv.writer(outFile)
		writer.writerow(['seq','count'])
		for seq in list(seqCount.keys()):
			writer.writerow([seq, seqCount[seq]])
	
	with open(options['files']['seqCountRange'],'w') as outFile:
		writer = csv.writer(outFile)
		writer.writerow(['seq','count'])
		for seq in list(seqCountRange.keys()):
			writer.writerow([seq, seqCountRange[seq]])

	# Create map of sequences 
	mapOfSeqs(options, seqCount, seqInit, refProt)



if __name__ == "__main__":
	main()
