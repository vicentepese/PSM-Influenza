### PSM-Narcolepsy
The 2009 Pandemrix (_PAN_) Influenza A pH1N1 vaccine has been linked to an increase of Narcolepsy type 1 (NT1) incidence in children. Such a phenomenon was not detected in similar vaccines like Arepandrix (_ARP_) or Foccetria (_FOC_). This repository was created to assess the relationship between Pandemrix and NT1,  through the study of MASSpect protein sequencing in the aforementioned vaccines to evaluate structural differences in the protein chain.

### Preparation 
The given data is filtered, taking only the instances of the _HA-_ proteins (_PAN_ and _ARP_) and _hemagglutinin_ (in _FOC_) defined in the 19<sup>th</sup> column of the given data. The name of the protein, the sequence, the initial position, and the vaccine name are taken. The sequence of peptides of reference is retrieved from [Entrez](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/)

### Settings and options
For the sake of maintaining an organized structure and simplicity across scripts, any option (such as directories, modifiable parameters, or other) is stored in `options.json`. The user must feel free to modified the file to adapt the repository locally, without affecting the functionality of the scripts here presented. Furthermore, certain options allow personalization, such as choosing the position range (which determines the positions analyzed referenced to the protein of reference). It is advised not to modify other files &ndash; unless absolutely necessary &ndash; to avoid any interruption of interaction between scripts. 

### PSM-Influenza
This section attempts to give a thorough description of the repository's content. 

### Outputs
This folder contains an instance of the outputs of the main functions. Since Github does not allow _HTML_ formatting, files must be downloaded and opened in an internet browser or _HTML_ reader. 

### Utils
The scripts contained in this folder are of general and/or common use across the main scripts of the repository. 
- `filterData.sh` filters and takes relevant information from the raw data. 
- `mergeDataMASS.py` merges the post-filtered data into a single file for easier manipulation and load. 
- `scrolling_template.html` is a _.html_ file that serves as an initial template for final Outputs. 
- `utils.py` contains common functions used in main scripts. For further information about the functions, it is advised to read the file. 

### v1
This folder contains older or previous versions of scripts, now unused. 

### Sequence mapping and counting

Each script is designed to perform a specific task, yet all share the same `options.json` as expressed previously. The functions displayed in this section only consider the _ARP_ and _PAN_ vaccines. 
- `mapSequences.py` maps every fragment or sequence of Aminoacids (AA) based on a range of positions defined in `options.json`, and counts the number of sequences (instances) of each vaccine. Outputs are written in a _.txt_ file.
- `mapSequencesHTML.py` is a more complex and functional version of `mapSequences.py`. It maps every fragment or sequence of AA based on a range of positions defined in `options.json` as well, and in addition, highlights PTMs using a unique color code and frames mutations (using the protein of reference previously retrieved). The script gives the percentage of the total of sequences for each sequence and vaccine. The output is in _.html_ format and must be opened with a browser.
- `mapPTM.py` counts the number of PTMs for each position and PTM type, and gives a percentage of the total of positions for each vaccine. The output is in _.html_ format and must be opened with a browser.
- `uniqueSeqsPTM.py` maps unique sequences and counts the number of PTMs per position and type, giving a percentage of the total of instances of the sequence. The script highlights the position of the PTM in <span color='red'>red</span> and present for each of the highlight position the type and proportion of PTMs following the sequence, and using a unique color code. 

### Jacob's presentation method
The scripts' outcomes presented in this section are inspired by the format used by Jacob L. et al. in their study ["Comparison of Pandemrix and Arepanrix, two pH1N1 AS03-adjuvanted vaccines differentially associated with narcolepsy development"](https://www.sciencedirect.com/science/article/pii/S0889159114005194?via%3Dihub). Scrips presented in this section study _PAN_ versus _ARP_ and _FOC_. 
- `jacobMethod.py` maps the post-translational modifications presented in _ARP_ and _FOC_ and tests their significance against _PAN_. A significance level of 0.05 was set and is reflected in _p-values_ in red. 
- `mutJaconMethod.py` maps the mutations presented in _ARP_, _FOC_ and _PAN_ in comparison to the protein of reference. A significance level of 0.05 was set. 
- `testGlycosylation.py` test the significant occurrence of at least one glycosylation between the positions 277 and 323. 
