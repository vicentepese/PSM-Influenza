### PSM-Narcolepsy
This repository was designed to assess the relation between Narcolepsy and the Influenza vaccine through the study of MASSpect protein sequences in two Influenza Virus Vaccine strains (_PAN_ and _APR_).

#### Preparation 
The given data is filtered, taking only the instances of the _HA-_ proteins (19<sup>th</sup> column of the given data). Only the name of the protein, the sequence, and the initial position are taken. The sequence of peptides of reference si retrieved from [Entrez](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/)

#### Sequence mapping and counting

Each script is desined to perform a specific task, yet all share the same `options.json` file for the sake of simplicity. 
- `mapSequences.py` maps every fragment or sequence of Aminoacids (AA) based on a range of positions defined in `options.json`, and counts the number of sequences (instances) of each vaccine. The outputs is written in a _.txt_ file.
- `mapSequencesHTML.py` is a more complex and functional version of `mapSequences.py`. It maps every fragment or sequence of AA based on a range of positions defined in `options.json` as well, and in addition, highlights PTMs using a unique color code and frames mutations (using the protein of reference previously retrieved). The script gives the percentage of the total of sequences for each sequence and vaccine. The output is in _.html_ format and must be opened with a browser.
- `mapPTM.py` counts the number of PTMs for each position and PTM type, and gives a percentage of the total of positions for each vaccine. The output is in _.html_ format and must be opened with a browser.
- `uniqueSeqsPTM.py` maps unique sequences and counts the number of PTMs per position and type, giving a percentage of the total of instances of the sequence. The script highlights the position of the PTM in <span color='red'>red</span> and present for each of the hightlight position the type and proportion of PTMs following the sequence, and using a unique color code. 


### mapSequences.py

First, the sequence of the filtered data (thus, only the *HA-* proteins) are mapped for each position between the range set in `options.json` ) &ndash; ergo, the sequence of Aminoacids (AA) is written for each position between the set position range. Only sequences that overlap at least one position with the defined position range are considered. Posttranslational Modifications (PTM), defined between squared brackets in the sequence (\[\]), have been removed of the sequence (not affecting the indexing of the positions). 

A map for each position within the range defined in `options.json` is created in order to find potential mutations. Then, each sequence is counted for both vaccines. Finally, sequences are mapped based on their initial position, and refered to the protein of reference. 

### mapSequencesHTML.py

The protein of reference is retrieved, and the pre-filtered MASSpect data corresponding to the _HA_ sequences of the vaccines is loaded. Sequences and Posttranslational Modifications (PTM) are counted. The vaccine sample is computed in order to further give a percentage of sequence ocurrence. Each sequence is subsequently and individually mapped; for each sequence, mutations are found, PTMs are located and counted. 

The output is a _.html_ file that highlights each PTM position and type using a color code, and frames mutations. The percentage of isntances for each sequence and vaccine is included. 

### mapPTM.py

Firstly, the protein of reference is retrieved and the pre-filtered data is loaded. Then, for each position in the protein of reference (and in accordance with the position range defined in `options.json`), PTMs are counted for each vaccine, giving a percentage of instances for each PTM and vaccine. 

The output is a _.html_ file that allows to visualize where PTMs are mostly taking place, and which is the frequency of occurence. 

### uniqueSeqsPTM.py

The script starts retrieving the protein of reference and loading the data. Then, PTMS are counted for each sequence, position and type. 

The output is a _.html_ file that highlights in  <span color='red'>red</span> the location of the PTMS, and for such locations includes the type of location found in each position, and the frequency of sequences that presented such PTM (for that sequence).







Finally, in order to visually assess the mutations, each sequence is chronologically sorted based on their initial positions, and compared to a reference protein (_refProt_ in `options.json`). The number of repetitions of the sequence for each strain is included at the end of the positions. Mutations have been highlighted
