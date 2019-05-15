---
layout: default-overview
title: Genome-Guided Transcriptome Assembly
exercises: 45
questions:
  - How to assemble my RNAseq genome-guided?
objectives:
  - Run the different tools
  - Understand each step of the assembly
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly
```

## Trimmomatic/Hisat2/Stringtie

### Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

```
cd $RNAseq_assembly_path
mkdir -p guided_assembly/trimmomatic
cd guided_assembly
module load trimmomatic/0.36
```

The following command line will perform the following:
	<br>• Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
	<br>• Remove leading low quality or N bases (below quality 3) (LEADING:3)
	<br>• Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
	<br>• Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
	<br>• Drop reads below the 36 bases long (MINLEN:36)

```
java -jar /sw/apps/bioinfo/trimmomatic/0.36/milou/trimmomatic-0.36.jar PE -threads 5 -phred33 $data/raw_computes/ERR305399_1.fastq.gz $data/raw_computes/ERR305399_2.fastq.gz trimmomatic/ERR305399.left_paired.fastq.gz trimmomatic/ERR305399.left_unpaired.fastq.gz trimmomatic/ERR305399.right_paired.fastq.gz trimmomatic/ERR305399.right_unpaired.fastq.gz ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```


### Hisat2

Once the reads have been trimmed, we use [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) to align the RNA-seq reads to a genome in order to identify exon-exon splice junctions.
HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) against a reference genome.

First you need to build an index of your genome

```
mkdir index

module load HISAT2/2.1.0
module load samtools/1.8
hisat2-build $data/genome/genome.fa index/genome_index
```

Then you can run Hisat2 :

**--phred33** Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the "Phred+33" encoding, which is used by the very latest Illumina pipelines (you checked it with the script fastq_guessMyFormat.pl)
<br>**--rna-strandness** For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF. (RF means fr-firststrand see [here](https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/rnaseq_library_types.md) for more explanation).
<br>**--novel-splicesite-outfile** In this mode, HISAT2 reports a list of splice sites in the file :
chromosome name tab genomic position of the flanking base on the left side of an intron tab genomic position of the flanking base on the right tab strand (+, -, and .) '.' indicates an unknown strand for non-canonical splice sites.

```
mkdir hisat2

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
module load StringTie/1.3.3

stringtie hisat2/accepted_hits.sorted.bam -o stringtie/transcripts.gtf
```

When done you can find your results in the directory ‘stringtie’. The file transcripts.gtf includes your assembled transcripts.

:bulb: **Tips**: You could now also visualise all this information using a genome browser, such as IGV. IGV requires a genome fasta file and any number of annotation files in GTF or GFF3 format (note that GFF3 formatted file tend to look a bit weird in IGV sometimes).

Transfer the gtf files to your computer using scp:

```
scp __YOURLOGIN__@rackham.uppmax.uu.se:/proj/g2019006/nobackup/__YOURLOGIN__/RNAseq_assembly/guided_assembly/stringtie/transcripts.gtf .
```

:question: Looking at your results, are you happy with the default values of Stringtie (which we used in this exercise) or is there something you would like to change?

:bulb: **Tips**: Maybe some statistics about the results would be nice, run the command

```
stringtie
```
To have all the parameters available.

<details>
<summary>:key: Click to see the solution .</summary>

<ul>If you want to have the gene abundance information for instance you should use the parameters -A </ul>
<ul>You can also use a reference annotation file if your genome has been annotated already and you want to use this annotation in your assembly -G </ul>
<ul>You can be more or less selective on the isoform abundance and keep really low abundant isoform or discard them -f </ul>
<ul>You can decide if you want to keep only reads with high coverage and set the minimum read coverage higher than the default parameter (2.5) -c </ul>
<br>There are many parameters to play with depending on your question.

Check the <a href="https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual/">Stringtie manual</a> for more information.

</details>

## Check the intron size of your genes (Optional)

From the gtf file, you can now know the size of introns in your genes.

:bulb: **Tips**: This can be an important step later when you are running the structural annotation and you need to write in the parameters what introns size the tool should be expecting (this parameter for Maker, augustus and genemark exists).

To do so, first you need to convert your gtf into a proper gff3 (You have done it in the Practical: Abinitio with augustus) and then run the script gff3_sp_manage_introns.pl

:bulb: **Tips**: Do
```
gxf_to_gff3.pl --help
gff3_sp_manage_introns.pl --help
```

<details>
<summary>:key: Click to see the solution .</summary>

<code>

gxf_to_gff3.pl -g stringtie/transcripts.gtf -o stringtie/transcript_stringtie.gff3
<br>gff3_sp_manage_introns.pl --gff stringtie/transcript_stringtie.gff3 -o introns_information

</code>
</details>

:question: What is the value you should choose?

<details>
<summary>:key: Click to see the solution .</summary>
You can choose 6722, 6500 if you think the introns size is overestimate or 7000 if you think it is possible to have bigger introns and you do not want to miss them.
</details>
