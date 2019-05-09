---
layout: default-overview
title: Assessing quality of you RNAseq assembly
exercises: 45
questions:
  - How do I check the quality of my assembly?
objectives:
  - Run busco on the de-novo and genome-guided assembly
  - Understand the output
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly
```

# Assessing the quality using busco

There are different ways of assessing the quality of your assembly, you will find some of them [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment).

We will run busco to check the the quality of the assembly.
[BUSCO](https://busco.ezlab.org/) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness (what we are going to do here). Genes that make up the BUSCO sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species in the chosen branch of tree of life.

For the trinity results :

```
cd $RNAseq_assembly_path

mkdir assembly_assessment

cd assembly_assessment

module load bioinfo-tools
module load BUSCO
source $BUSCO_SETUP

run_BUSCO.py -i $data/RNAseq/trinity/Trinity.fasta -o busco_trinity -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5
```

Busco will take 30 min to run so you can check the results in $data/RNAseq/busco_trinity


For the guided assembly results

You need first to extract the transcript sequences from the gtf transcript file :

```
ln -s $RNAseq_assembly_path/guided_assembly/stringtie/transcripts.gtf .

gff3_sp_extract_sequences.pl --cdna -g transcripts.gtf -f $data/genome/genome.fa -o transcripts_stringtie.fa
```
Then you can run busco again :

```
run_BUSCO.py -i stringtie/transcripts_stringtie.fa -o busco_stringtie -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5
```

Compare the two busco, what do you think happened for stringtie?


# What's next?

Now you are ready either to annotate your RNAseq or you can use then to do the genome annotation.

For the de-novo assembly you can use the Trinity.fasta file obtained.
For the genome-guided assembly you can either use the Stringtie results transcripts.gtf but you will often need to reformat it into a gff file.
