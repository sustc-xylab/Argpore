# What ARGpore Does
ARGpore is designed to predict resistome of 2D.fasta/1D.fasta generated from NanoPore sequencing. ARGpore can help users to identify ARG hits on nanopore reads and simultaneously find the carrier population of the identified ARGs by searching against phylogenetic marker genes on nanopore reads

NOTE: You may still use ARGpore for your metagenomic-assembled contig/scaffold, but you would lose bouche of ARGs and taxa due to annotation limitation.


*Citation:*

If you use ARGpore in your nanopore dataset analysis please cite:

Xia, Yu, An-Dong Li, Yu Deng, Xiao-Tao Jiang, Li-Guan Li, and Tong Zhang. MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. Frontiers in Microbiology 2017


## Pre-requisites for installation

	lastal	there is a version of last-744 will be downloaded with argpore in this git 
		please make sure lastal and lastdb in last-744/scr are execuatable and in your PATH
		before run argpore.sh 
	ruby
	R

## Installation 

git clone Argpore to your pc

	git clone https://github.com/sustc-xylab/Argpore.git

	cd Argpore

## Download MetaPhlAn marker gene database at below link:

https://1drv.ms/u/s!AuLWV06xGWsWhfIGqUedoF3gX-CAzA

or this link:

https://pan.baidu.com/s/1K6eFfgZRWJSyBi8dRG1AfA

## Build marker gene index

	lastdb -Q 0 markers.lastindex markers.fasta -P 10

## Run argpore simply like this:

	bash argpore.sh -f test.fa 

## Detailed usage:

*Input:*

input of ARGpore is simply your 2D.fasta/1D.fasta

Resistome prediction Algorithm:

STEP1: Your input fasta will be searched against nt-version SARG database (v1.0)

STEP2: Valid alignment (with > 60% similarity over > 90% alignment length by default) will be kept for further filtering of overlap regions

STEP3: If two hit regions on the same read overlaped for > 50% alignment length, only the one with longest ARG hit will be kept,

*Taxa Algorithm:*

STEP1: Your input fasta will be searched against clade specific marker gene database (MetaPhlan 2)

STEP2: Valid alignment (with > 60% similarity over > 90% alignment length by default) will be kept for taxa annotation

STEP3: Only the best hit with highest bitscore is kept to determine phylogenetic affiliation

*Output:*

input_taxa.tab:nanopore reads with valid match to clade specific marker gene

input_arg.tab:nanopore reads with valid match to SARG database

input_arg.w.taxa.tab: ARG-containing nanopore reads with taxa annotated
