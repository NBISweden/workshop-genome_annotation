---
layout: default-overview
title: De-novo Transcriptome Assembly
exercises: 20
questions:
  - How to De-novo assemble my RNAseq?
  - What should I look in my assembly to go forward
objectives:
  - run Trinity
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly
```

## Trinity

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) assemblies can be used as complementary evidence, particularly when trying to polish a gene build with Pasa. Before you start, check how big the raw read data is that you wish to assemble to avoid unreasonably long run times.

```
cd $RNAseq_assembly_path

mkdir trinity

cd trinity

module load trinity/2.4.0
module load samtools/1.9

Trinity --seqType fq --max_memory 32G --left $data/raw_computes/ERR305399_1.fastq.gz --right $data/raw_computes/ERR305399_2.fastq.gz --CPU 5 --output trinity --SS_lib_type RF
```

Trinity takes a long time to run (several hours), you can stop the program when you start it and have a look at the results, look in $data/RNAseq/trinity the output is Trinity.fasta

:bulb: **Tips**: Using Trinity in genome annotation

Some advantages :
- Really easy to run.
- Many manual an tutorial available on internet eg : [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki).
- Output can be annotated directly with [Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki).
- Widely use in the community.
- One can concatenate several read files/libraries or run them separately.

Some disadvantage :
- Can take several hours to run.
- Can take a lot of resources to run and create a lot of data.
- Input data are often really big files :
As a "rule" in our team, if read files are >10Gb compressed, one should consider normalizing the data. If you have such big file, it is recommended to normalize them prior to assembly to avoid very long run times.
- Results tend to be noisy with many shorts transcripts, might need to be filtered.
