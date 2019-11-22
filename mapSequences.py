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

class Format:
	underline = '\033[0m]'
	end = '\033[4m]'

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

				seqMap[str(pos)].append(AA)
		
		elif end_pos >= pos_range[0] and end_pos <= pos_range[1]:

			# Map each position 
			for pos, AA in zip(range(init_pos, (end_pos + 1)), AAseq):

				seqMap[str(pos)].append(AA)
			
	return seqMap

def countSequences(options, AFLSA):

	# For each sequence 
	pos_range = options['pos_range']

	# TODO: Count sequences
	seqCount = defaultdict(int)
	seqCountRange = defaultdict(int)
	seqInit = defaultdict(int)
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

			# Count 
			seqCount[AAseq] += 1

			# Count in range 
			seqCountRange[AAseq[0:(pos_range[1]-init_pos)]] += 1

			# Keep initial position 
			seqInit[AAseq] = init_pos
		
		elif end_pos >= pos_range[0] and end_pos <= pos_range[1]:

			# Count 
			seqCount[AAseq] += 1

			# Count in range 
			seqCountRange[AAseq[pos_range[0] - init_pos:]] += 1

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
	
	# Convert to string index 
	if len(pos_mut) > 0:
		pos_mut_idx = [pos-init_pos for pos in pos_mut]
	else:
		pos_mut_idx = []

	return pos_mut_idx


def mapOfSeqs(options, seqCount, seqInit, refProt):

	# Create array of positions: min initial pos to max ending pos
	init_pos = min(seqInit[seq] for seq in list(seqInit.keys()))
	ending_pos = max(seqInit[seq] + len(seq) for seq in list(seqInit.keys()))
	indx_array = np.arange(init_pos,ending_pos+1)

	# Indexed protein string 
	pos_range = options['pos_range']
	refProt_string = [refProt[pos] for pos in list(refProt.keys()) if (pos >= pos_range[0]) & (pos <= pos_range[1])]
	refProt_string = ' '*(pos_range[0] - init_pos) + ''.join(refProt_string) + ' '*(ending_pos - pos_range[1])

	# Create string for each sequence 
	seqString = defaultdict(str)
	seqMut = defaultdict(str)
	for seq, seq_init_pos in list(seqInit.items()):
		seqString[seq] = ' '*(np.absolute(init_pos-seq_init_pos)) + seq + '(' + str(seqCount[seq]) + ')'

		# Find mutations, get index for the string
		pos_mut_idx = findMutations(refProt, seq, seq_init_pos)
		
		# Keep position for string
		if len(pos_mut_idx) > 0:
			pos_mut_idx.insert(0,0)
			mut_idx_array = [pos_mut_idx[i] - pos_mut_idx[i-1] for i in range(1,len(pos_mut_idx))]
			seq_mut_string = ' '*(np.absolute(init_pos-seq_init_pos))
			for mut_idx in mut_idx_array:
				seq_mut_string = seq_mut_string + ' '*(mut_idx) + '_'
			seqMut[seq] = seq_mut_string
		else:
			seqMut[seq] = ''

	# Order by initial position
	seqInit = OrderedDict([(k, seqInit[k]) for k in sorted(seqInit, key=seqInit.get)])

	# Write 
	with open('seqMap.txt','w') as outFile:
		outFile.writelines(str(idx)[0] for idx in indx_array)
		outFile.write('\n')
		outFile.writelines(str(idx)[1] for idx in indx_array)
		outFile.write('\n')
		outFile.writelines(str(idx)[2] for idx in indx_array)
		outFile.write('\n')
		outFile.write('\n')
		outFile.write(refProt_string + '\n')
		outFile.write('\n')
		for seq in list(seqInit.keys()):
			outFile.write(seqMut[seq] + '\n')
			outFile.write(seqString[seq] + '\n')


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