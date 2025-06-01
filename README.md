# bat_project

![alt text](https://github.com/Malaevleo/bat_project/blob/main/git_title.jpg "Title Card")


*!!! Some parts of the pipeline are not present in the pipeline.ipynb. It is explicitely stated when and how these steps should be done !!!*

# Data preparation part

**MMSEQS2 homologs search** \
First things first, our analysis begins with finding out which genes to analyze. After extracting human genes from databases (in our case we used OpenGenes and The Matrisome Project) we run MMSEQS2 tool to find out homologs in bats. In order to do this task in a convenient way we merge all bat proteomes into one big fasta file and assign tags to each protein based on species of origin. This will allow us to later easily identify protein of origin and create a file with best matches between target and query protein. 

Run on your server:
```
mmseqs/bin/mmseqs easy-search human_genes.fasta bat_proteomes_merged.fasta results.m8 -c 0.2 tmp --min-seq-id 0.3 
```
After obtaining results of the homologs search we need to create file where best matches for each human protein and each bat species are stored. 

**Biopython paired fasta files generation**

In order to continue analysis we have to create separate files for each protein consisting of nucleotide sequences from all of the bats. Additionally, we should have files with amino acid sequences as they will be useful for structural modeling, for example.

**Alignment using MACSE and conversion using PAL2NAL**

MACSE is an essential part as it creates proper alignments while accounting for codon structure of genes. To use it run the script ```run_macse.sh```. \
PAL2NAL is used to convert MACSE alignments into proper codon alignments suitable for the $\frac{dN}{dS}$ analysis. To use it run the script ```run_pal2nal.sh```.

**Trees pruning and labeling**
Trees should be labeled for each protein so that HYPHY understands phylogeny when analysing each gene. Foreground and background branches should also be labeled. In our case foreground branches were the ones with long-living bats.

# dN/dS HYPHY

Scheme for the pipeline:
![alt text](https://github.com/Malaevleo/bat_project/blob/main/pipeline%20adj.jpg "Pipeline Scheme")

First step of the dN/dS analysis is BUSTED-E as it is the most robust tool for gene-level positive selection identification yet [1]. Initially we run analysis on the Foreground (**run_bustede1.sh**), then, for genes experiencing positive selection, we run BUSTED-E on the Background branches (**run_bustede2.sh**) in order to check whether this event also happens within the gene on the other branches of the tree. 

Resulting logs are analyzed in the pipeline.ipynb with applied FDR correction. Genes undergoing positive selection can later be used for the other steps of the pipeline. You can copy only significant ones into other folder (and, if you want, delete existing one that will become obsolete after this step) using **copy_sig.sh**.

After finding out which genes experience positive selection we can double-check our prediction using BUSTED-PH. To run BUSTED-PH use **run_bustedph.sh**. If you have little time or severely limited computational resources, you can skip this part.

To find out at which branches positive selection happened we can employ branch-level model such as aBSREL. To run aBSREL use **run_absrel.sh**.

In order to obtain 'root-to-tip' $\frac{dN}{dS}$ values suitable for the PGLS analysis we have to run FitMG94 with ```--type lineage``` option. To run FitMG94 use **run_mg94.sh**. 

# PGLS

Run it using BAT1K.R file. Additional data prepocessing steps might be required.

# TODO
1. Docker version
2. Single pipeline file for server usage
3. Tests folder
