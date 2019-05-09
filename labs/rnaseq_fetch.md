---
layout: default-overview
title: Preparing RNAseq data for assembly
exercises: 45
questions:
  - How do I check the quality of my RNAseq data?
objectives:
  - Run fastqc
  - Understand the output
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/sw/courses/annotation/2019/data
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly
mkdir -p $RNAseq_assembly_path
```

# Introduction

This exercise is meant to get you acquainted with the type of data you would normally encounter in an annotation project. We will for all exercises use data for the fruit fly, Drosophila melanogaster, as that is one of the currently best annotated organisms and there is plenty of high quality data available.

You can create the folders where you want but I would suggest a folder organisation, if you do not follow this organisation remember to put the correct path to your data

```
cd $RNAseq_assembly_path

```

## Assembling transcripts based on RNA-seq data

Rna-seq data is in general very useful in annotation projects as the data usually comes from the actual organism you are studying and thus avoids the danger of introducing errors caused by differences in gene structure between your study organism and other species.

Important remarks to remember before starting working with RNA-seq:

- Check if RNAseq are paired or not. Last generation of sequenced short reads (since 2013) are almost all paired. Anyway, it is important to check that information, which will be useful for the tools used in the next steps.

- Check if RNAseq are stranded. Indeed this information will be useful for the tools used in the next steps. (In general way we recommend to use stranded RNAseq to avoid transcript fusion during the transcript assembly process. That gives more reliable results. )

- Left / L / forward / 1 are identical meaning. It is the same for Right / R /Reverse / 2


### Checking encoding version and fastq quality score format

To check the technology used to sequences the RNAseq and get some extra information we have to use fastqc tool.

```
module load bioinfo-tools
module load FastQC/0.11.5

mkdir fastqc_reports

fastqc $data/raw_computes/ERR305399_1.fastq.gz -o fastqc_reports/
```
scp the html file resulting of fastqc, what kind of result do you have?

```
scp __YOURLOGIN__@rackham.uppmax.uu.se:/proj/g2019006/nobackup/__YOURLOGIN__/RNAseq_assembly/fastqc_reports/YOURFILE .
```
Checking the fastq quality score format :

As we will be using the scripts libraries available in the git gaas you need first to export the libraries (if you were logged off):

Then :
```
fastq_guessMyFormat.pl -i $data/raw_computes/ERR305399_1.fastq.gz

```

In the normal mode, it differentiates between Sanger/Illumina1.8+ and Solexa/Illumina1.3+/Illumina1.5+.
In the advanced mode, it will try to pinpoint exactly which scoring system is used.

More test can be made and should be made on RNA-seq data before doing the assembly, we have not time to do all of them during this course. have a look [here](https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools)

ADD JACQUES/MASTER STUDENTS SCRIPT TO CHECK RNA?
