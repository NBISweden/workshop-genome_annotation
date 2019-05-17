---
layout: default-overview
title: Preparing RNAseq data for assembly
exercises: 30
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
export data=/proj/g2019006/nobackup/$USER/data
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

RNA-seq data is in general very useful in annotation projects as the data usually comes from the actual organism you are studying and thus avoids the danger of introducing errors caused by differences in gene structure between your study organism and other species.

Important remarks to remember before starting working with RNA-seq:

- Check if RNAseq are paired or not. Last generation of sequenced short reads (since 2013) are almost all paired. Anyway, it is important to check that information, which will be useful for the tools used in the next steps.

- Check if RNAseq are stranded. Indeed this information will be useful for the tools used in the next steps. (In general way we recommend to use stranded RNAseq to avoid transcript fusion during the transcript assembly process. That gives more reliable results. )

- Left / L / forward / 1 are identical meaning. It is the same for Right / R /Reverse / 2


### Checking encoding version and fastq quality score format

To check the technology used to sequences the RNAseq and get some extra information we have to use fastqc tool.

```
module load FastQC/0.11.5

mkdir fastqc_reports

fastqc $data/raw_computes/ERR305399_1.fastq.gz -o fastqc_reports/
```
Copy the html file resulting of fastqc in your local machine

```
scp __YOURLOGIN__@rackham.uppmax.uu.se:/proj/g2019006/nobackup/__YOURLOGIN__/RNAseq_assembly/fastqc_reports/YOURFILE .
```
:question: what kind of results do you have?

<details>
<summary>:key: Click to see the solution .</summary>
<br>
Fastqc reports give you different statistics about your RNAseq data before assembly.
Next to each categories, there is color code to tell you when data is good (green), bad or missing information (red) and questionable results (orange).  

<br>You can find more details about the results <a href="https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/">here</a>.
</details>

<br>
Checking the fastq quality score format :

```
fastq_guessMyFormat.pl -i $data/raw_computes/ERR305399_1.fastq.gz
```

In the normal mode, it differentiates between Sanger/Illumina1.8+ and Solexa/Illumina1.3+/Illumina1.5+.
In the advanced mode, it will try to pinpoint exactly which scoring system is used.

More test can be made and should be made on RNA-seq data before doing the assembly, we do not have time to do all of them during this course. Have a look [here](https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools)

### Checking the library type (read orientations)

The information regarding library type can be very useful for reads to be assembled into a transcriptome or mapped to a reference assembly. This is because the library type can help to discern where in the transcriptome shorter ambiguous reads belong by using the read’s relative orientation and from which strand it was sequenced. Unfortunately, this information regarding the library type used is not included in sequencing output files and may be lost before the assembly of the data.

Here a resume of the different library types:  

 <img align="center" src="https://raw.githubusercontent.com/NBISweden/GUESSmyLT/master/library_types.jpg"  />

In order to guess de-novo the library type we will use for this excercise: [GUESSmyLT](https://github.com/NBISweden/GUESSmyLT).  

First link the data and load the necessary modules:  
```bash
cd $RNAseq_assembly_path
ln -s $data/raw_computes/ERR305399_1.fastq.gz
ln -s $data/raw_computes/ERR305399_2.fastq.gz
ln -s $data/genome/genome.fa

module load trinity/2.8.2
module load BUSCO/3.0.2b
source $BUSCO_SETUP
module load bowtie2/2.3.4.3
module load samtools/1.2
module load python3/3.6.0
module load snakemake/5.4.5
```

Then install [GUESSmyLT](https://github.com/NBISweden/GUESSmyLT):   
```
git clone https://github.com/NBISweden/GUESSmyLT.git
cd GUESSmyLT
python3 setup.py install --user
```

Check that the tool is working fine:  
```bash
~/.local/bin/GUESSmyLT -h
```

Let's know launch it on our data:  
```bash
cd $RNAseq_assembly_path
 ~/.local/bin/GUESSmyLT --reads ERR305399_1.fastq.gz ERR305399_2.fastq.gz --reference genome.fa --threads 10 --subsample 100000 --output result
```

:question:
<ol><li>Do you recognize the pipeline framework run by GUESSmyLT? </li>
<li>What library type has been used for the ERR305399 sample? </li>

<details>
<summary>:key: Click to see the solution .</summary>
<br>
<ol>
<li>GUESSmyLT run <strong>Snakemake</strong>.</li>
<li>The data is <strong>fr-unstranded</strong> (according to the <code>first strand</code>) or in other term <strong>RF unstranded</strong> (according to <code>mRNA</code>).</li>
</ol>
</details>
