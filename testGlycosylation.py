import sys
import json
import re
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from scipy import stats
from utils import getBindingCore, importBindData,\
	 importData, reference_retreive, getBindingCore


def getGlycoAmid(options, data):

    # Initialize 
    vaccSample = defaultdict(int)
    glycoCount = defaultdict(int)
    nonglycoCount = defaultdict(int)
    nonglycoSample = defaultdict(int)
    nonDeamidCount = defaultdict(int)
    nonDeamidSample = defaultdict(int)
    nonGlDa = defaultdict(int)
    nonGlDaSample = defaultdict(int)
    glycoRange = options['glycoRange']

    # Count the number of glycosylation
    for seq in data:

        # Initialize sequence
        AAseq = seq[1][2:-2]
        AAnonPTM = re.sub('\[.+?\]', '', AAseq)
        initPos = int(seq[2])
        endPos = initPos + len(AAnonPTM) -1
        pos_range = options['pos_range']
        NPos = 277 - initPos


        # If segement overlaps with defined range
        if not(initPos >= pos_range[1]) and not(endPos <= pos_range[0]): 
            
            # Count sample 
            vaccSample[seq[3]] +=1

            # Count if glycosylated 
            if '+203.079' in AAseq:
                glycoCount[seq[3]] += 1
            
            ########### SCENARIO B ###########
            # If sequence contains both 277N and glycosylation sites
            if not(initPos >= glycoRange[1] ) and not(endPos <= glycoRange[0]) or 277 in range(initPos, endPos+1):
                
                # Find indexes and instances of PTMs
                PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
                PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

                # Initialize variables
                idx_cumm = 0
                flag = 0

                # For each instance, if there is a glycosylation within the range, raise a flag
                for instance, idx in zip(PTM_instances, PTM_idx):
                    if instance == '+203.079' and idx.start() -1 -idx_cumm + initPos in range(glycoRange[0], glycoRange[1]+1):
                        flag = 1
                    idx_cumm += len(instance) + 2
                
                # If there is a deamidatio in 277N, raise a flag
                if not(AAseq[NPos+2:NPos+8] == '+0.984'):
                    flag = 1

                # If flag was not raised, then the sequence did not contain neither deamidation nor glycosylation. Count. 
                if flag == 0:
                    nonGlDa[seq[3]] += 1
                
                # Count the sample of sequences that contained both 277N and glycosylation sites
                nonGlDaSample[seq[3]] +=1

            ########### SCENARIO A ###########

            # If the sequence contains glycsylation site
            if not(initPos >= glycoRange[1] ) and not(endPos <= glycoRange[0]):

                # Find indexes and instances of PTMs
                PTM_idx = re.finditer('\[(.*?)\]', AAseq, re.DOTALL)
                PTM_instances = re.findall('\[(.*?)\]', AAseq, re.DOTALL)

                # Initialize variables
                idx_cumm = 0
                flag = 0

                # For each instance, check glycosylation and raise flag
                for instance, idx in zip(PTM_instances, PTM_idx):
                    if instance == '+203.079' and idx.start() -1 -idx_cumm + initPos in range(286, 288+1):
                        flag = 1
                    idx_cumm += len(instance) + 2
                
                # If flag not raised (no glyco.) count
                if flag == 0:
                    nonglycoCount[seq[3]] += 1  

                # Count sample
                nonglycoSample[seq[3]] += 1                          
                    
            # Count non-amidated 277N
            if 277 in range(initPos, endPos+1):
                if not(AAseq[NPos+2:NPos+8] == '+0.984'):
                    nonDeamidCount[seq[3]] += 1
                nonDeamidSample[seq[3]] += 1
            

    return vaccSample, glycoCount, nonglycoCount, nonglycoSample, nonDeamidCount, nonDeamidSample, nonGlDa, nonGlDaSample

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

    # Get glycosylation count, and non-glycosylated sites and non-amidated 277N 
    vaccSample, glycoCount, nonglycoCount, nonglycoSample \
        , nonDeamidCount, nonDeamidSample, nonGlDa, nonGlDaSample = getGlycoAmid(options, data)

    # Compute fisher test
    fisherResult = fisherTest(glycoCount, vaccSample)

    # Print 
    print('\n')
    print('Glycosylation Fisher\'s test result:')
    for vacc in fisherResult:
        print(vacc + ' {:.2%} ({})'.format(glycoCount[vacc]/vaccSample[vacc], vaccSample[vacc]) + \
            ' vs. PAN {:.2%} ({}): pvalue = {:.2}, oddsratio = {:.2}'.format( \
            glycoCount['PAN']/vaccSample['PAN'], vaccSample['PAN'], \
            fisherResult[vacc]['pvalue'], fisherResult[vacc]['oddsratio']))

    # Print probabilities of deamidation or glycosylation as separate entities which can 
    # apper in different sequences. This method does not allow to compute a Fisher's test, but 
    print('\n')
    print('SCENARIO A: Considering 277N and glycosylation sites as different entities that can occur in separate segments')
    print('Probability of non-deamydated and non-glycosylated in ARP:' + \
         '{:.2%}'.format((nonglycoCount['ARP']/nonglycoSample['ARP']) * (nonDeamidCount['ARP']/nonDeamidSample['ARP'])))
    print('Probability of non-deamydated and non-glycosylated in FOC:' + \
         '{:.2%}'.format((nonglycoCount['FOC']/nonglycoSample['FOC']) * (nonDeamidCount['FOC']/nonDeamidSample['FOC'])))     
    print('Probability of non-deamydated and non-glycosylated in PAN:' + \
         '{:.2%}'.format((nonglycoCount['PAN']/nonglycoSample['PAN']) * (nonDeamidCount['PAN']/nonDeamidSample['PAN'])))

    # Print probabilities of sequences that include both 277N and glycosylation sites. 
    # This allows to compute Fisher's test PAN vs ARP and FOC.  
    print('\n')
    print('SCENARIO B: Considering only sequence that included both 277N and glycosylation sites:')
    print('Probability of non-deamydated and non-glycosylated in ARP:' + \
          '{:.2%}'.format((nonGlDa['ARP']/nonGlDaSample['ARP'])))
    print('Probability of non-deamydated and non-glycosylated in FOC:' + \
          '{:.2%}'.format((nonGlDa['FOC']/nonGlDaSample['FOC'])))
    print('Probability of non-deamydated and non-glycosylated in PAN:' + \
          '{:.2%}'.format((nonGlDa['PAN']/nonGlDaSample['PAN'])))

     # Fisher test on SCENARIO B:
    nonGlDaFisher = fisherTest(nonGlDa, nonGlDaSample)

    # Print 
    print('\n')
    for vacc in nonGlDaFisher:
        print('Fisher\'s test result for SCENARIO B:')
        print(vacc + ' {:.2%} ({})'.format(nonGlDa[vacc]/nonGlDaSample[vacc], nonGlDaSample[vacc]) + \
            ' vs. PAN {:.2%} ({}): pvalue = {:.2}, oddsratio = {:.2}'.format( \
            nonGlDa['PAN']/nonGlDaSample['PAN'], nonGlDaSample['PAN'], \
            nonGlDaFisher[vacc]['pvalue'], nonGlDaFisher[vacc]['oddsratio']))



        

if __name__ == "__main__":
    main()
