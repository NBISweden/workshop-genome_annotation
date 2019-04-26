---
layout: default
title:  'Exercise RNAseq assembly'
---

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](uppmax_login).

# Assessing quality of you RNAseq assembly

There are different ways of assessing the quality of your assembly, you will find some of them [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment).

We will run busco to check the the quality of the assembly.
[BUSCO](https://busco.ezlab.org/) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness (what we are going to do here). Genes that make up the BUSCO sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species in the chosen branch of tree of life.

For the trinity results :

```
cd ~/annotation_course/RNAseq_assembly

mkdir busco

cd busco

module load BUSCO/3.0.2b

/sw/apps/bioinfo/BUSCO/3.0.2b/rackham/bin/run_BUSCO.py -i ~/annotation_course/RNAseq_assembly/trinity/Trinity.fasta -o busco_trinity -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5
```

Busco will take 30 to run so you can check the results in ~/annotation_course/????/RNAseq/busco_trinity


For the guided assembly results

You need first to extract the transcript sequences from the gtf transcript file :

```
module load BioPerl
export PERL5LIB=$PERL5LIB:~/RNAseq_assembly_annotation/GAAS/annotation/
~/RNAseq_assembly_annotation/GAAS/annotation/Tools/bin/gff3_sp_extract_sequences.pl --cdna -g transcripts.gtf -f ~/RNAseq_assembly_annotation/assembly_annotation/chromosome/chr4.fa -o transcripts_stringtie.fa

```
Then you can run busco again :

```
/sw/apps/bioinfo/BUSCO/3.0.2b/rackham/bin/run_BUSCO.py -i stringtie/transcripts_stringtie.fa -o busco_stringtie -l $BUSCO_LINEAGE_SETS/arthropoda_odb9 -m tran -c 5

```

Compare the two busco, what do you think happened for stringtie?


# What's next?

Now you are ready either to annotate your RNAseq or you can use then to do the genome annotation.

For the de-novo assembly you can use the Trinity.fasta file obtained.
For the genome-guided assembly you can either use the Stringtie results transcripts.gtf but you will often need to reformat it into a gff file.
