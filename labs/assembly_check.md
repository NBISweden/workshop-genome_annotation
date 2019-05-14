---
layout: default-overview
title: Assembly Check
exercises: 30
questions:
  - How do I check the quality of my assembly regarding gene contents?
  - What should I look in my assembly to go forward
objectives:
  - run busco and understand the output
  - have a look at the fasta and its statistics
---

# Introduction

Before starting an annotation project, we need to carefully inspect the assembly to identify potential problems before running expensive computes.
You can look at i) the fragmentation (N50, N90, how many short contigs); ii) the sanity of the fasta file (presence of Ns, presence of ambiguous nucleotides, presence of lowercase nucleotides, single line sequences vs multiline sequences); iii) completeness using BUSCO; iv) presence of organelles; v) others (GC content, how distant the investigated species is from the others annotated species available).
The two next exercices will perform some of these checks.

# Prerequisites  
For this exercise you need to be logged in to Uppmax.

Setup the environment:
```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export assembly_check_path=/proj/g2019006/nobackup/$USER/assembly_check
mkdir -p $assembly_check_path
```

# 1 Checking the gene space of your assembly

[BUSCO](https://busco.ezlab.org/) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness. Genes that make up the BUSCO sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species.

***Note:*** In a real-world scenario, this step should come first and foremost. Indeed, if the result is under your expectation you might be required to enhance your assembly before to go further.

:mortar_board: **_Exercise 1_ - BUSCO -:**

You will run BUSCO on the genome assembly.

First create a busco folder where you work:
```bash
cd $assembly_check_path
mkdir busco
cd busco
```

The [BUSCO website](http://busco.ezlab.org) provides a list of datasets containing the cores genes expected in the different branches of the tree of life. To know in which part/branch of the tree of life is originated your species, you can look at the [NCBI taxonomy website](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7227) (Lineage line).
You have to chose the most appropriate BUSCO dataset for your species. You may download the latest BUSCO dataset from the [busco website](http://busco.ezlab.org) or use the one already present at uppmax (`ll $BUSCO_LINEAGE_SETS/`).

Now you are ready to launch BUSCO to check the completness of your assembly (genome.fa file).

```
module load bioinfo-tools
module load BUSCO
source $BUSCO_SETUP
run_BUSCO.py -i $data/genome/genome.fa -o genome_dmel_busco -m genome -c 10 -l $BUSCO_LINEAGE_SETS/arthropoda_odb9
```

While BUSCO is running, you may start the exercise 2 (to do so you will need to open another terminal).
When done, check the short\_summary\_genome\_dmel\_busco file in the output folder.

:question: How many core genes have been searched in your assembly ? How many are reported as complete? Does this sound reasonable?  
:bulb: **Tips**: the "genome" is here in fact only the chromosome 4 that corresponds to less than 1% of the real size of the genome.

# 2 Various Checks of your Assembly

:mortar_board: **_Exercise 2_ :**
Launching the following script will provide you some useful information.

```
cd $assembly_check_path
fasta_statisticsAndPlot.pl -f $data/genome/genome.fa -o fasta_checked
```

:question: Is your genome very fragmented (number of sequences)? Do you have high GC content ? Do you have lowercase nucleotides ? Do you have N at sequence extremities?

If you don't see any peculiarities, you can then decide to go forward and start to perform your first wonderful annotation.
