# HotWASp 
A tool to boost GWAS results through the utilization of networks. 

## Installation
Install HotWASp by cloning into the repository:
```
git clone https://github.com/aryakaul/HotWASp.git
cd HotWasp
python3 scripts/hotwasp.py --help
```

## How to run
`HotWASp` has two required arguments: `-n, --network` and `-g, --gwas`. 

#### Input Network
The `-n` argument takes in a full path to a network for HotWASp to use. This network should take the following file format "GeneName\tGeneName\n". As a result, the file should have as many lines as edges present. An example of a file is produced below:

```
head -n 4 examplenetwork
```

```
GeneA   GeneB
GeneA   GeneC
GeneD   GeneC
GeneD   GeneA
```

Note that we do not take in any information on edge-weight. If your initial network is a weighted network, we recommend thresholding the edge weights at a specific cutoff, *c*, and only including them in the file if the edge weight is > *c*. 


#### Input GWAS
The `-g` argument takes in a full path to a file containing GWAS results for HotWASp to use. These results should take the following file format: "GeneName\tP-value\n". Where every gene tested is accompanied with its respective adjusted p-value. An example of a file is produced below:

```
head -n 4 examplegwas
```

```
GeneA 0.0003
GeneB 0.8
GeneC 0.02
GeneD 1e-10
```


