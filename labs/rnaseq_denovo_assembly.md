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
<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](uppmax_login).

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/sw/courses/annotation/2019/data
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly

```

## Trinity

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) assemblies can be used as complementary evidence, particularly when trying to polish a gene build with Pasa. Before you start, check how big the raw read data is that you wish to assemble to avoid unreasonably long run times.


PATH

```
cd $RNAseq_assembly_path

mkdir trinity

cd trinity

module load bioinfo-tools
module load trinity/2.4.0
module load samtools

Trinity --seqType fq --max_memory 32G --left $data/raw_computes/ERR305399_1.fastq.gz --right $data/raw_computes/ERR305399_2.fastq.gz --CPU 5 --output trinity --SS_lib_type RF
```

Trinity takes a long time to run (like hours), you can stop the program when you start it and have a look at the results, look in ~/RNAseq_assembly_annotation/assembly_annotation/RNAseq/trinity the output is Trinity.fasta
