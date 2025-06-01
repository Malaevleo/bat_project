# bat_project

![alt text](https://github.com/Malaevleo/bat_project/blob/main/pipeline%20adj.jpg "Pipeline Scheme")

Blastp results for human matrisome against genome of Myotis Myotis (resultsmyofiltered.txt): 

```
makeblastdb -in myotis_protein.faa -dbtype prot -out myotis_myotis_proteome 
blastp -query matrisome/refseq_sequences.fasta -db myotis_myotis_proteome -out resultsmyofiltered.txt -outfmt 6 -num_threads 6 -max_target_seqs 3
```

# Data preparation part

**MMSEQS2 homologs search** \
First things first, our analysis begins with finding out which genes to analyze. After extracting human genes from databases (in our case we used OpenGenes and The Matrisome Project) we run MMSEQS2 tool to find out homologs in bats. In order to do this task in a convenient way we merge all bat proteomes into one big fasta file and assign tags to each protein based on species of origin. This will allow us to later easily identify protein of origin and create a file with best matches between target and query protein. 

```
mmseqs/bin/mmseqs easy-search human_genes.fasta bat_proteomes_merged.fasta results.m8 -c 0.2 tmp --min-seq-id 0.3 
```
After obtaining results of the homologs search we need to create file where best matches for each human protein and each bat species are stored. 

**Biopython paired fasta files generation**

In order to continue analysis we have to create separate files for each protein consisting of nucleotide sequences from all of the bats. Additionally, we should have files with amino acid sequences as they will be useful for structural modeling, for example.

**Alignment using MACSE and conversion using PAL2NAL**

MACSE is an essential part as it creates proper alignments while accounting for codon structure of genes. \
PAL2NAL is used to convert MACSE alignments into proper codon alignments suitable for the $\frac{dN}{dS}$ analysis.

**Trees pruning and labeling**
Trees should be labeled for each protein so that HYPHY understands phylogeny when analysing each gene. Foreground and background branches should also be labeled. In our case foreground branches were the ones with long-living bats.

# dN/dS HYPHY

