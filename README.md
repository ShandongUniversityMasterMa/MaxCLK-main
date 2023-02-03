# MaxCLK: discovery of cancer driver genes via maximal clique and information entropy of modules

### This is the original repository for the MaxCLK paper. 

**Requirements**

Ubantu 21.04

```
<fstream>
<vector>
<string>
<iostream>
<stdlib.h>
<algorithm>
<cmath>
<eigen3/Eigen/Dense>
<math.h>
<ctime>
<map>
```

## **Input**

### Protein Protein Interaction (PPI) network

There are two kinds of files: the Protein Protein Interaction (PPI) network edges file and the nodes to index mapping file, including three different PPI networks. All the PPI network files can be download from http://compbio-research.cs.brown.edu/pancancer/hotnet2/

1.The PPI network file

For example:
The file is located at PPR_data/PPR_hi.m. The row (or column) of the element corresponds to the index of the gene.

```
0.402912139892578,4.31300804848433e-06,3.79393941329909e-06, ...
4.31299486081116e-06,0.445848286151886,1.23754898595507e-05, ...
1.26464590266551e-06,4.12517420045333e-06,0.40184274315834, ...
...
```

2. The nodes to index file mapping

For example:
The file is located at PPR_data/hi_index_file.txt
```
index Gene_Name
1 A1BG
2 A1CF
3 A2BP1
...
```

### Mutation Data

1.Patients and their mutant genes

For example:
The file is located at mutation_data/LAML

```
TCGA-AB-2970 BCLAF1 BCOR FLT3 PLEKHH1 RUNX1 WT1
TCGA-AB-2971 E2F8 MAPK1 SAP130 TET2 WEE1
TCGA-AB-2972 C10orf118 DROSHA FANCI LARP4B MORC3 NPM1 PTPN11 SRSF6 STAG2 TUFT1 ZC3H18 ZNF43
....
```

### reference driver genes

This file is located at ref_genes/NCGgenes.txt


## **Run**

There are fourteen bash scripts corresponding to the different mutation datasets, including three pan-cancer datasets and eleven sub-datasets of a single cancer type.

For example:

```
cd ./runSH/

bash ./LAML.sh
```

## **Outputs**

For cancer type LAML, MaxCLK will output the resulting module information, including the gene composition in each module and the IE (p_value of IE) of the module.

```
./out_module/LAML/hi/e0.1_l0_u3/CSKmodule_3.txt

Coverage_in_minimun_set_cover Coverage Cov Ex CovEx IE p_value_IE gene1 gene2 gene3

87 87 0.457895 1 0.457895 0.842039 0.00968367 NPM1 RUNX1 TP53
21 29 0.152632 1 0.152632 0.97415 0.00774693 CEBPA KIT KRAS
16 32 0.168421 1 0.168421 0.84266 0.00968367 KDM6A TET2 WT1
...
```

```
./out_module/LAML/hi/e0.1_l0_u3/CSKmodule_3_genes.txt

gene1 gene2 gene3

NPM1 RUNX1 TP53
CEBPA KIT KRAS
KDM6A TET2 WT1
...
```

The consensus result will be stored at consensus/e0.1_LAML/subnet_A_s.txt. (s=0,1,...,7)
