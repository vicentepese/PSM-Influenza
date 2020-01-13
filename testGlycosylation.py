import sys
import json
import re
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from scipy import stats
from utils import getBindingCore, importBindData,\
	 importData, reference_retreive, getBindingCore


def getGlyco(options, data):

    # Initialize 
    glycoCount = defaultdict(int)
    vaccSample = defaultdict(int)

    # Count the number of glycosylation
    for seq in data:

        # Initialize sequence
        AAseq = seq[1][2:-2]
        initPos = int(seq[2])
        endPos = initPos + len(AAseq)
        pos_range = options['pos_range']

        # If segement within defined range
        if not(initPos > pos_range[1]) and not(endPos < pos_range[0]): 

            # Count
            if '+203.079' in AAseq:
                glycoCount[seq[3]] += 1

            # Count sample 
            vaccSample[seq[3]] +=1

    return glycoCount, vaccSample

def fisherTest(glycoCount, vaccSample):

    # Initialize
    fisherResult = defaultdict(lambda : defaultdict(int))
    
    # For each vaccine compute Fisher test
    for vacc in ['ARP','FOC']:

        # Create arrays 
        glycoPos = [glycoCount[vacc], glycoCount['PAN']]
        glycoNeg = [vaccSample[vacc] - glycoCount[vacc], \
            vaccSample['PAN'] - glycoCount['PAN']]

        # Compute test 
        oddsratio, pvalue = stats.fisher_exact([glycoPos, glycoNeg])
        fisherResult[vacc]['oddsratio'] = oddsratio
        fisherResult[vacc]['pvalue'] = pvalue
    
    return fisherResult

def main():

    # Read options
    with open('options.json') as inFile:
        options = json.load(inFile)
    
    # Set pos range to analyze segement of interest
    options['pos_range'] = [277, 323]

    # Import data 
    data = importData(options)

    # Count the number of glycosylation within the segment
    glycoCount, vaccsample = getGlyco(options, data)

    # Compute fisher test
    fisherResult = fisherTest(glycoCount, vaccsample)

    # Print 
    for vacc in fisherResult:
        print(vacc + ' {:.2%} ({})'.format(glycoCount[vacc]/vaccsample[vacc], vaccsample[vacc]) + \
            ' vs. PAN {:.2%} ({}): pvalue = {:.2}, oddsratio = {:.2}'.format( \
            glycoCount['PAN']/vaccsample['PAN'], vaccsample['PAN'], \
            fisherResult[vacc]['pvalue'], fisherResult[vacc]['oddsratio']) + \
                 '\n')

if __name__ == "__main__":
    main()
