---
layout: default-overview
title: Genome-Guided Transcriptome Assembly
exercises: 65
questions:
  - How to assemble my RNAseq genome-guided?
objectives:
  - Run the different tools
  - Understand each step of the assembly
---

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](uppmax_login).

## Trimmomatic/Hisat2/Stringtie

### Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

```
cd ~/annotation_course/RNAseq_assembly
mkdir -p guided_assembly/trimmomatic
module load bioinfo-tools
module load trimmomatic/0.36
```

The following command line will perform the following:
	• Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
	• Remove leading low quality or N bases (below quality 3) (LEADING:3)
	• Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
	• Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
	• Drop reads below the 36 bases long (MINLEN:36)

UPDATE PATH

```
java -jar /sw/apps/bioinfo/trimmomatic/0.36/milou/trimmomatic-0.36.jar PE -threads 5 -phred33 ~/annotation_course/data/raw_computes/ERR305399_1.fastq.gz ~/annotation_course/data/raw_computes/ERR305399_2.fastq.gz trimmomatic/ERR305399.left_paired.fastq.gz trimmomatic/ERR305399.left_unpaired.fastq.gz trimmomatic/ERR305399.right_paired.fastq.gz trimmomatic/ERR305399.right_unpaired.fastq.gz ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.36/milou/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

```


### Hisat2

Once the reads have been trimmed, we use [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) to align the RNA-seq reads to a genome in order to identify exon-exon splice junctions.
HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) against a reference genome.

First you need to build an index of your chromosome 4

```
mkdir index

module load HISAT2
module load samtools/1.8
hisat2-build ~/annotation_course/data/genome/genome.fa index/genome.fa_index
```

Then you can run Hisat2 :
--phred33 Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the "Phred+33" encoding, which is used by the very latest Illumina pipelines (you checked it with the script fastq_guessMyFormat.pl)
--rna-strandness <string> For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF. (RF means fr-firststrand see [here](https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/rnaseq_library_types.md) for more explanation).
--novel-splicesite-outfile <path> In this mode, HISAT2 reports a list of splice sites in the file :
chromosome name <tab> genomic position of the flanking base on the left side of an intron <tab> genomic position of the flanking base on the right <tab> strand (+, -, and .) '.' indicates an unknown strand for non-canonical splice sites.

```
hisat2 --phred33 --rna-strandness RF --novel-splicesite-outfile hisat2/splicesite.txt -S hisat2/accepted_hits.sam -p 5 -x index/genome_index -1 trimmomatic/ERR305399.left_paired.fastq.gz -2 trimmomatic/ERR305399.right_paired.fastq.gz
```

Finally you need to change the sam into bam file and to sort it in order for stringtie to use the bam file to assemble the read into transcripts.

```
samtools view -bS -o hisat2/accepted_hits.bam hisat2/accepted_hits.sam

samtools sort -o hisat2/accepted_hits.sorted.bam hisat2/accepted_hits.bam
```


### Stringtie

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus.
You can add as input an annotation from gtf/gff3 file to calculate TPM and FPKM values.


```
module load bioinfo-tools
module load StringTie

stringtie hisat2/accepted_hits.sorted.bam -o stringtie/transcripts.gtf
```

When done you can find your results in the directory ‘outdir’. The file transcripts.gtf includes your assembled transcripts.

You could now also visualise all this information using a genome browser, such as IGV. IGV requires a genome fasta file and any number of annotation files in GTF or GFF3 format (note that GFF3 formatted file tend to look a bit weird in IGV sometimes).

Transfer the gtf files to your computer using scp:

```
scp __YOURLOGIN__@rackham.uppmax.uu.se:~/annotation_course/RNAseq_assembly/stringtie/transcripts.gtf .
```

Looking at your results, are you happy with the default values of Stringtie (which we used in this exercise) or is there something you would like to change?

CHANGE GTF TO GFF?
