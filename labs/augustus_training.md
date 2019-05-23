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
export structural_annotation_path=/proj/g2019006/nobackup/$USER/structural_annotation
export augustus_training_path=/proj/g2019006/nobackup/$USER/structural_annotation/augustus_training
mkdir -p $augustus_training_path
```

# Introduction

From the maker run evidence based, we can train our ab-initio predictors and then use them for the second run of annotation.
You will need a set of genomic sequences with gene structures (sequence coordinates of starts and ends of exons and genes) and the most important part is to select the right set of genes.
In many cases, or as a first step towards modeling complete genes, it is sufficient to have only the coding parts of the gene structure (CDS).
We will only train augustus today as it is one the best ab-initio predictor and one of the hardest to train.
Maker also support SNAP (Works good, easy to train, not as good as others ab-initio especially on longer intron genomes), GeneMark (Self training, no hints, buggy, not good for fragmented genomes or long introns), FGENESH (Works great, costs money even for training) and now EVM.


## Training Augustus

Move into the folder where we will train Augustus
```
cd $augustus_training_path
```

You will need to symlink the evidence-based annotation (the gff annotation file from the first run of maker) and the genome fasta sequence.
```
ln -s $structural_annotation_path/maker/maker_evidence/maker.gff maker_evidence.gff
ln -s $data/genome/genome.fa
```
## Compile a set of training and test genes

First step is to select only the coding genes from the maker.gff file and remove all tRNA ( :bulb: **Tips**: in this case there no tRNA but it is important to remove them)
```
gff3_sp_splitByLevel2Feature.pl -g maker_evidence.gff -o maker_results_noAbinitio_clean
ln -s maker_results_noAbinitio_clean/mrna.gff
```
In this folder you will need to create different folders
```
mkdir filter  
mkdir protein  
mkdir nonredundant  
mkdir blast_recursive  
mkdir gff2genbank  
```
Next step, we need to filter the best genes we will use for the training, we need complete genes and a AED score under 0.3 (:bulb: **Tips**:those are our default parameters you can change them if you want to be more selective with an AED under 0.1, you can also set a distance between genes if you know your genes are further appart).

```
maker_select_models_by_AED_score.pl -f mrna.gff -v 0.3 -t "<=" -o filter/codingGeneFeatures.filter.gff
```

We also need to select the longest ORF (we want to use the complete longest genes to avoid creating incorrect and too short genes when we use the abinitio profile )
```
gff3_sp_keep_longest_isoform.pl -f filter/codingGeneFeatures.filter.gff -o filter/codingGeneFeatures.filter.longest_cds.gff
```
/!\\ If you receive this error "Nothing to do... this file doesn't contain any isoform !", it means that you do not have any isoforms. Thus, you can directly use filtered annotation (based on AED) for following step.
```
cp filter/codingGeneFeatures.filter.gff filter/codingGeneFeatures.filter.longest_cds.gff
```

There are different ways of proceeding after the first selection and we are using "the approached of spliced alignments of protein sequences of the same or a very closely related species" against the assembled genomic sequence.
In order to do so, we translate our coding genes into proteins, format the protein fasta file to be able to run a recursive blast and then select the best ones.
Indeed, each sequence can contain one or more genes; the genes can be on either strand. However, the genes must not overlap, and only one transcript per gene is allowed.
```
gff3_sp_extract_sequences.pl -g filter/codingGeneFeatures.filter.longest_cds.gff -f genome.fa -o protein/codingGeneFeatures.filter.longest_cds.proteins.fa

module load blast/2.7.1+   

makeblastdb -in protein/codingGeneFeatures.filter.longest_cds.proteins.fa -dbtype prot  

blastp -query protein/codingGeneFeatures.filter.longest_cds.proteins.fa -db protein/codingGeneFeatures.filter.longest_cds.proteins.fa -outfmt 6 -out blast_recursive/codingGeneFeatures.filter.longest_cds.proteins.fa.blast_recursive

gff3_sp_filter_by_mrnaBlastValue_bioperl.pl --gff filter/codingGeneFeatures.filter.longest_cds.gff --blast blast_recursive/codingGeneFeatures.filter.longest_cds.proteins.fa.blast_recursive --outfile nonredundant/codingGeneFeatures.nr.gff

```
Sequences need to be converted to a simple genbank format.
```
module load augustus/3.2.3

gff2gbSmallDNA.pl nonredundant/codingGeneFeatures.nr.gff $data/genome/genome.fa 500 gff2genbank/codingGeneFeatures.nr.gbk
```
In order for the test accuracy to be statistically meaningful the test set should also be large enough (100-200 genes).
You should split the set of gene structures randomly.
```
module load cufflinks/2.2.1
randomSplit.pl gff2genbank/codingGeneFeatures.nr.gbk 100
```
:question:What happened? how can you solve it? what might be the consequences of it?

<details>
<summary>:key: Click to see the solution .</summary>
There are not 100 genes in the file, because we are using only the chr4 of drosophila.
The training will probably not be good!
</details>

## Train Augustus

Now that you have created a set of gene to train augustus, let's train it!

```
module load BUSCO/3.0.2b

source $BUSCO_SETUP

new_species.pl --species=dmel_$USER

etraining --species=dmel_$USER gff2genbank/codingGeneFeatures.nr.gbk.train

augustus --species=dmel_$USER gff2genbank/codingGeneFeatures.nr.gbk.test | tee run.log 
```
- Look at the accuracy report, what does it mean? why? see [Training Augustus](http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html)
