---
layout: default-overview
title: Assessing quality of your RNAseq assembly
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

### For the trinity results :

```
cd $RNAseq_assembly_path

mkdir assembly_assessment

cd assembly_assessment

module load BUSCO/3.0.2b
source $BUSCO_SETUP

run_BUSCO.py -i $data/RNAseq/trinity/Trinity.fasta -o busco_trinity -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5
```

Busco will take 30 min to run so you can check the results in $data/RNAseq/busco_trinity


### For the guided assembly results:

You need first to extract the transcript sequences from the gtf transcript file :

```
ln -s $data/genome/genome.fa
gff3_sp_extract_sequences.pl --cdna -g $RNAseq_assembly_path/guided_assembly/stringtie/transcripts.gtf -f genome.fa -o $RNAseq_assembly_path/guided_assembly/stringtie/transcripts_stringtie.fa
```
Then you can run busco again :

```
run_BUSCO.py -i $RNAseq_assembly_path/guided_assembly/stringtie/transcripts_stringtie.fa -o busco_stringtie -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5
```

:question:Compare the two busco, what do you think happened for stringtie?

<details>
<summary>:key: Click to see the solution .</summary>
We only used the chromosome 4 of the Drosophila as genome to do the assembly with stringtie while Trinity is not mapped to any chromosome and so contain all the transcripts for the complete genome. BUSCO compares a set of genes of a complete genome and not only a part of it.
It makes no sense to use BUSCO on only 1 chromosome of a genome! :)

:bulb:Also if you notice there are many duplicates in the BUSCO results, In this case it is due to the fact that all isoforms have been kept so each isoform is consider as 1 gene. You need to select one of them (like the longest for instance with the script gff3_sp_keep_longest_isoform.pl or the one you prefer) and you will have a more accurate results of BUSCO.
</details>

# What's next?

Now you are ready use the results of your De-novo assembly and guided assembly to do the genome annotation.

For the de-novo assembly you can use the Trinity.fasta file obtained.
For the genome-guided assembly you can either use the Stringtie results transcripts.gtf but you will often need to reformat it into a gff file.
If you have not done it please do :
```
gxf_to_gff3.pl -g $RNAseq_assembly_path/guided_assembly/stringtie/transcripts.gtf -o $RNAseq_assembly_path/guided_assembly/stringtie/transcript_stringtie.gff3
```

You are now ready to use the genome-guided assembly for your annotation.
