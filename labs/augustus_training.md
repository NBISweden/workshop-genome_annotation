---
layout: default-overview
title: Training ab-initio predictor
exercises: 60
questions:
  - What do I need to train Augustus?
  - What are the step to train Augustus?
objectives:
  - Going through the steps
  - Understanding the different part of the training (it is complex so take your time!)
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export augustus_training_path=/proj/g2019006/nobackup/$USER/structural_annotation/augustus_training
mkdir -p $augustus_training_path
```

# Introduction

From the maker run evidence based, we can train our ab-initio predictors and then use them for the second run of annotation.
You will need a set of genomic sequences with gene structures (sequence coordinates of starts and ends of exons and genes) and the most important part is selected the right set of genes.
In many cases, or as a first step towards modeling complete genes, it is sufficient to have only the coding parts of the gene structure (CDS).
We will only train augustus today as it is one the best ab-initio predictor and one of the hardest to train.
Maker also support SNAP (Works good, easy to train, not as good as others ab-initio especially on longer intron genomes), GeneMark (Self training, no hints, buggy, not good for fragmented genomes or long introns), FGENESH (Works great, costs money even for training) and now EVM.


## Training Augustus

cd into the folder and then load all modules that we will need to train Augustus
```
cd $augustus_training_path

module load bioinfo-tools   
module load cufflinks/2.2.1
```

You will need to symlink all data you will need such as the gff files from the first run of maker and the chromosome 4 fasta sequence.
```
ln -s $structural_annotation_path/maker/maker_evidence/maker.gff maker_no_abinitio.gff
```
## Compile a set of training and test genes

First step is to select only the coding genes from the maker.gff file and remove all tRNA ( :bulb: **Tips**: in this case there no tRNA but it is important to remove them)
```
gff3_sp_splitByLevel2Feature.pl -g maker_no_abinitio.gff -o maker_results_noAbinitio_clean

cd maker_results_noAbinitio_clean
```
In this folder you will need to create different folders
```
mkdir filter  
mkdir protein  
mkdir nonredundant  
mkdir blast_recursive  
mkdir gff2genbank  
```
Next step, we need to filter the best genes we will use for the training, we need complete genes and a AED score under 0.3 (:bulb: **Tips**:those are our default parameters you can change them if you want to be more selective with an AED under 0.1, you can alsoknow your gene are further appart or if you want to be more selective with an AED under 0.1).

```
maker_select_models_by_AED_score.pl -f mrna.gff -v 0.3 -t "<=" -o filter/codingGeneFeatures.filter.gff
```

We also need to select the longest ORF (we want to use the complete longest genes to avoid create incorrect and too short genes when we use the abinitio profile )
```
gff3_sp_keep_longest_isoform.pl -f filter/codingGeneFeatures.filter.gff -o filter/codingGeneFeatures.filter.longest_cds.gff
```
There are different ways of proceeding after the first selection and we are using "the approached of spliced alignments of protein sequences of the same or a very closely related species" against the assembled genomic sequence.
In order to do so, we translate our coding genes into proteins, format the protein fasta file to be able to run a recursive blast and then select the best ones.
Indeed, each sequence can contain one or more genes; the genes can be on either strand. However, the genes must not overlap, and only one transcript per gene is allowed.
```
gff3_sp_extract_sequences.pl -g filter/codingGeneFeatures.filter.longest_cds.gff3 -f $data/genome/genome.fa -o protein/codingGeneFeatures.filter.longest_cds.proteins.fa

module load blast  

makeblastdb -in protein/codingGeneFeatures.filter.longest_cds.proteins.fa -dbtype prot  

blastp -query protein/codingGeneFeatures.filter.longest_cds.proteins.fa -db protein/codingGeneFeatures.filter.longest_cds.proteins.fa -outfmt 6 -out blast_recursive/codingGeneFeatures.filter.longest_cds.proteins.fa.blast_recursive

gff3_sp_filter_by_mrnaBlastValue_bioperl.pl --gff codingGeneFeatures.filter.longest_cds.gff --blast blast_recursive/codingGeneFeatures.filter.longest_cds.proteins.fa.blast_recursive --outfile nonredundant/codingGeneFeatures.nr.gff

```
Sequences need to be converted in a simple genbank format.
```
gff2gbSmallDNA.pl nonredundant/codingGeneFeatures.nr.gff 4.fa 500 gff2genbank/codingGeneFeatures.nr.gbk
```
In order for the test accuracy to be statistically meaningful the test set should also be large enough (100-200 genes).
You should split the set of gene structures randomly.
```
randomSplit.pl gff2genbank/codingGeneFeatures.nr.gbk 100
```
- What happened? how can you solve it? what might be the consequences of it?


## Train Augustus

Now that you have created a set of gene to train augustus, let's train it!

Augustus need a set of parameters that are provided :

please use the path where you copied augustus_path in the Busco exercise yesterday.
```
module load augustus/3.2.3

new_species.pl --AUGUSTUS_CONFIG_PATH=augustus_path --species=dmel_login

AUGUSTUS_CONFIG_PATH=augustus_path

etraining --species=dmel_login gff2genbank/codingGeneFeatures.nr.gbk.train 

augustus --species=dmel_login gff2genbank/codingGeneFeatures.nr.gbk.test | tee run.log 
```
- Look at the accuracy report, what does it mean? why? see [Training Augustus](http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html)
