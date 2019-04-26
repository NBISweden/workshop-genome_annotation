---
layout: default
title:  'Exercise RNAseq assembly'
---

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](uppmax_login).

## De-novo Transcriptome Assembly

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) assemblies can be used as complementary evidence, particularly when trying to polish a gene build with Pasa. Before you start, check how big the raw read data is that you wish to assemble to avoid unreasonably long run times.


PATH

```
cd ~/annotation_course/RNAseq_assembly

mkdir trinity

cd trinity

module load bioinfo-tools
module load trinity/2.4.0
module load samtools

Trinity --seqType fq --max_memory 32G --left ~/RNAseq_assembly_annotation/assembly_annotation/raw_computes/ERR305399_1.fastq.gz --right ~/RNAseq_assembly_annotation/assembly_annotation/raw_computes/ERR305399_2.fastq.gz --CPU 5 --output trinity --SS_lib_type RF
```

Trinity takes a long time to run (like hours), you can stop the program when you start it and have a look at the results, look in ~/RNAseq_assembly_annotation/assembly_annotation/RNAseq/trinity the output is Trinity.fasta


NECESSARY ???  not convinced :

In order to compare the output of stringtie and the output of trinity we need to map the trinity transcript to the chr4 of Drosophila.

We'll use the GMAP software to align the Trinity transcripts to our reference genome. Trinity contains a utility that facilitates running GMAP, which first builds an index for the target genome followed by running the gmap aligner:

```
module load gmap-gsnap

mkdir gmap

/sw/apps/bioinfo/trinity/2.4.0/rackham/util/misc/process_GMAP_alignments_gff3_chimeras_ok.pl --genome ~/RNAseq_assembly_annotation/assembly_annotation/chromosome/chr4.fa --transcripts ~/RNAseq_assembly_annotation/assembly_annotation/RNAseq/trinity/Trinity.fasta > gmap/transcript_trinity.gff
```
