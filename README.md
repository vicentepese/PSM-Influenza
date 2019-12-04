### PSM-Narcolepsy
This repository was designed to assess the relation between Narcolepsy and the Influenza vaccine through the study of protein sequences in two Influenza Virus Vaccine strains and protein sequences in the body.

#### Preparation 
The given data is filtered, taking only the instances of the _HA-_ proteins (19<sup>th<\sup> column of the given data). Only the name of the protein, the sequence, and the initial position are taken. The sequence of peptides of reference si retrieved from _[Entrez]_(https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/)

#### Sequence mapping and counting

Each script is desined to perform a specific task, yet all share the same `options.json` file for the sake of simplicity. 
- `mapSequences.py` maps every fragment or sequence of Aminoacids (AA) based on a range of positions defined in `options.json`, and counts the number of sequences (instances) of each vaccine. The outputs is written in a _.txt_ file.
- `mapSequencesHTML.py` is a more complex and functional version of `mapSequences.py`. It maps every fragment or sequence of AA based on a range of positions defined in `options.json` as well, and in addition, highlights PTMs using a unique color code and frames mutations (using the protein of reference previously retrieved). The script gives the percentage of the total of sequences for each sequence and vaccine. The output is in _.html_ format and must be opened with a browser.
- `mapPTM.py` counts the number of PTMs for each position and PTM type, and gives a percentage of the total of positions for each vaccine. The output is in _.html_ format and must be opened with a browser.
- `uniqueSeqsPTM.py` maps unique sequences and counts the number of PTMs per position and type, giving a percentage of the total of instances of the sequence. The script highlights the position of the PTM in <span color='red'>red</span> and present for each of the hightlight position the type and proportion of PTMs following the sequence, and using a unique color code. 


### mapSequences.py

First, the sequence of the filtered data (thus, only the *HA-* proteins) are mapped for each position between the range set in `options.json` ) "---" ergo, the sequence of Aminoacids (AA) is written for each position between the set position range. Only sequences that overlap at least one position with the defined position range are considered. Posttranslational Modifications (PTM), defined between squared brackets in the sequence (\[\]), have been removed of the sequence (not affecting the indexing of the positions). 

A map for each position within the range defined in `options.json` is created in order to find potential mutations, and positions more pr


Finally, in order to visually assess the mutations, each sequence is chronologically sorted based on their initial positions, and compared to a reference protein (_refProt_ in `options.json`). The number of repetitions of the sequence for each strain is included at the end of the positions. Mutations have been highlighted

#### Commentary
All the settings of the repository must be set in `options.json`. The storage of the data is of free-will, yet paths must be defined in the `options.json` file.

`mapSequences.py` outputs a _.txt_ file, whereas `mapSequencesHTML.py` outputs an _HTML_ file. `options.json` allows to modify the way the mutations are highlighted under _"MarkdownHigh"_.

### Results 
Results are available in the folder `./Outputs/`, both in _.txt_ format and _.html_. In order to properly visualize the results, please open the proper in a text editor, and the latter in a browser. 
