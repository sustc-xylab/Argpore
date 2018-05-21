# Argpore
# installation
	lastal, there is one version of last-744 came with argpore 
	ruby

# example syntex
# run argpore simply like this:
argpore.sh -f test.fa 

What ARGpore Does

ARGpore can help users to identify ARG hits on nanopore reads and simultaneously find the carrier population of the identified ARGs by searching against phylogenetic marker genes on nanopore reads

*Citation:*
if you use ARGpore in your nanopore dataset analysis please cite:
Xia, Yu, An-Dong Li, Yu Deng, Xiao-Tao Jiang, Li-Guan Li, and Tong Zhang. MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. Frontiers in Microbiology 2017

*ARGpore docs:*

ARGpore is designed to predict resistome of 2D.fasta/1D.fasta generated from NanoPore sequencing.

NOTE: You may still use ARGpore for your metagenomic-assembled contig/scaffold, but you would lose bouche of ARGs and taxa due to annotation limitation.

*Input:*

input of ARGpore is simply your 2D.fasta/1D.fasta

Resistome prediction Algorithm:

STEP1: Your input fasta will be searched against nt-version SARG database (v1.0)

STEP2: Valid alignment (with > 80% similarity over > 70% alignment length by default) will be kept for further filtering of overlap regions

STEP3: If two hit regions on the same read overlaped for > 50% alignment length, only the one with longest ARG hit will be kept,

*Taxa Algorithm:*

STEP1: Your input fasta will be searched against clade specific marker gene database (MetaPhlan 2)

STEP2: Valid alignment (with > 80% similarity over > 70% alignment length by default) will be kept for taxa annotation

STEP3: Only the best hit with highest bitscore is kept to determine phylogenetic affiliation

*Output:*

input_runtime_taxa.tab:nanopore reads with valid match to clade specific marker gene

input_runtime_arg.tab:nanopore reads with valid match to SARG database

input_runtime_arg.w.taxa.tab: ARG-containing nanopore reads with taxa annotated

# the defult similarity cutoff （0.8） and alignment length cutoff （0.7） is set to ensure > 50% exact base match to ARG reference sequence or marker gene as 0.8*0.7*0.9 (the average accurary of 1D nanopore reads） = 0.504 

