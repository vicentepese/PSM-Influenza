import subprocess
import re
import json
import csv 

def importData(options):

	# Import data  as list of lists
	data = list() 
	with open(options['files']['mergedFiltMASS'], 'r') as inFile:
		reader = csv.reader(inFile)
		for row in reader:
			data.append(row)
	
	return data

def netMHCIIpan(options, AAseq):
	
	# Retrieve allele and run netMHCIIpan
	allele = options['netMHCIIpan-3.2']['allele']
	subprocess.run(['bash', './run_netMHCIIpan.sh', allele])

	# Open file 
	with open(options['files']['tmpOut.out'],'r') as inFile:
		output = list()
		for row in inFile:
			output.append(row)
	output = output[13].split()
	affinity = output[8]
	bindingLevel = output[-1][2:]
	bindingCore = output[5]

	return affinity, bindingLevel, bindingCore

def getBindingCore(options, data):

	for seq in data:

		# Create temporal fasta file 
		with open(options['files']['tmpSeq.fasta'],'w') as outFile:
			outFile.write('>seq1\n')
			outFile.write(seq)
		
		# Run netMHCIIpan with sequence 
		affinity, bindingLevel, bindingCore = netMHCIIpan(options, seq)
		print("Binding core prediction task successfully completed")
		coreIdx = re.search(bindingCore, seq).span()

	
	return affinity, bindingLevel, bindingCore, coreIdx

def main():

	# Read options 
	with open('options.json', 'r') as inFile:
		options = json.load(inFile)

	# Import  data
	data = importData(options)

	# Get binding core of every sequence 
	getBindingCore(options, data)

	

if __name__ == "__main__":
	main()