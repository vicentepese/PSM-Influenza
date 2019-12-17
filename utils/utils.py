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


def div0(n, d):
	return n / d if d and n else 0


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

def getRandomColor(options, **kwargs):

	# Create color dictionnary
	if "PTM_count" in kwargs:
		PTM_count = kwargs['PTM_count']
		r = lambda: random.randint(75, 200)
		color = {PTM: ['<span style=\"background: '+'#%02X%02X%02X; font-weight: bold' %
			(r(), r(), r()) + '\">', "</span>"] for PTM in list(PTM_count.keys())}
	color['ARP'] = ['<span style=\"color: #800000; font-weight: bold; \">', '</span>']
	color['PAN'] = ['<span style=\"color: #000782; font-weight: bold; \">', '</span>']
	color['strongBinder'] = ['<span style=\"background: #0F9D58; font-weight: bold; \">', '</span>']
	color['weakBinder'] = ['<span style=\"background: #F4B400; font-weight: bold; \">', '</span>']
	color['red'] = ['<span style=\"color: red; font-weight: bold; \">', '</span>']
	color['mut'] = ['<span style=\"font-weight: bold\">', '</span>']
	
	return color
