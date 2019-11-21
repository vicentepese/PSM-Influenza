### PSM-Narcolepsy
This repository was designed to assess the relation between Narcolepsy and the Influenza vaccine through the study of protein sequences in the Influenza Virus Vaccine and protein sequences in the body.

#### Preparation 
The given data is filtered, taking only the instances of the _HA-_ proteins (19<sup>th<\sup> column of the given data). Only the name of the protein, the sequence, and the initial position are taken

#### Sequence mapping and counting
First, the sequence of the filtered data (thus, only the *HA-* proteins) are mapped for each position between the range set in `options.json` ) "---" ergo, the sequence of Aminoacids (AA) is written for each position between the set position range. Only sequences that overlap at least one position with the defined position range are considered. Furthermore, Posttranslational Modifications (PTM), defined between squared brackets in the sequence (\[\]) have been removed of the sequence (not affecting the indexing of the positions). 

Subsequently, the number of repeated sequences overlapping to some extent with the proposed range of position is counted. Similarly to mapping the sequences, only the AA sequences that overlap at least in one position with the defined range are considered, and PTM was removed. Two outputs are defined: `seqCount` and `seqCountRange`. The propper corresponds to the number of sequences overlapping at least one position with the defined range, considering the positions outside of the range. The latter considers the sequences that overlap at least one position with the range, but only the AAs within the range were considered (thus does not include any AA outside the defined range).

#### Commentary
All the settings of the repository must be set in `options.json`. The storage of the data is of free-will, yet paths must be defined in the `options.json` file. 