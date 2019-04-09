# Assembling transcripts based on RNA-seq data

Rna-seq data is in general very useful in annotation projects as the data usually comes from the actual organism you are studying and thus avoids the danger of introducing errors caused by differences in gene structure between your study organism and other species.

Important remarks to remember before starting working with RNA-seq:
- Check if RNAseq are paired or not. Last generation of sequenced short reads (since 2013) are almost all paired. Anyway, it is important to check that information, which will be useful for the tools used in the next steps.
- Check if RNAseq are stranded. Indeed this information will be useful for the tools used in the next steps. (In general way we recommend to use stranded RNAseq to avoid transcript fusion during the transcript assembly process. That gives more reliable results. )
- Left / L / forward / 1 are identical meaning. It is the same for Right / R /Reverse / 2

First create a dedicated folder to work in:
```
cd ~/annotation_course/practical2
mkdir RNAseq
cd RNAseq
```

## 1. Genome guided transcriptome assembly: 

### Checking encoding version and fastq quality score format

To check the technology used to sequences the RNAseq and get some extra information we have to use fastqc tool.

```
mkdir fastqc
cd fastqc
mkdir fastqc_reports
fastqc ~/annotation_course/data/RNAseq/fastq/ERR305399.left.fastq.gz -o fastqc_reports/
```

Transfer the html file resulting of fastqc to your computer using scp in another terminal:   
```
scp -i ~/.ssh/azure_rsa student@__IP__:/home/student/annotation_course/practical1/augustus/augustus_drosophila.gff .
```
Open it. What kind of result do you have?

Checking the fastq quality score format

```
fastq_guessMyFormat.pl -i ~/annotation_course/data/RNAseq/fastq/ERR305399.left.fastq.gz
```

In the normal mode, it differentiates between Sanger/Illumina1.8+ and Solexa/Illumina1.3+/Illumina1.5+.
In the advanced mode, it will try to pinpoint exactly which scoring system is used.

More test can be made and should be made on RNA-seq data before doing the assembly, we have not time to do all of them during this course. have a look [here](https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools)

### Trimmomatic (trimming reads)

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

The following command line will perform the following:  
	• Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)  
	• Remove leading low quality or N bases (below quality 3) (LEADING:3)  
	• Remove trailing low quality or N bases (below quality 3) (TRAILING:3)  
	• Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)  
	• Drop reads below the 36 bases long (MINLEN:36)  

```
cd ~/annotation_course/practical2/RNAseq
mkdir trimmomatic
cd trimmomatic

java -jar trimmomatic-0.32.jar PE -threads 8 ~/annotation_course/data/RNAseq/fastq/ERR305399.left.fastq.gz ~/annotation_course/data/RNAseq/fastq/ERR305399.right.fastq.gz ERR305399.left_paired.fastq.gz ERR305399.left_unpaired.fastq.gz ERR305399.right_paired.fastq.gz ERR305399.right_unpaired.fastq.gz ILLUMINACLIP:trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Tophat (splice-aware mapping reads to genome)

Once the reads have been trimmed, we use [tophat](https://ccb.jhu.edu/software/tophat/index.shtml) to align the RNA-seq reads to a genome in order to identify exon-exon splice junctions. It is built on the ultrafast short read mapping program [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml).

```
cd ~/annotation_course/practical2/RNAseq
mkdir tophat
cd tophat

tophat --library-type=fr-firststrand ~/annotation_course/data/genome/4.fa ../trimmomatic/ERR305399.left_paired.fastq.gz ../trimmomatic/ERR305399.right_paired.fastq.gz -p 8
```

This step will take a really long time so you can use the bam file located here ~/annotation_course/data/RNAseq/tophat/accepted_hits.bam

### Stringtie (Assembling reads into transcripts)

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. Its input can include not only the alignments of raw reads used by other transcript assemblers, but also alignments longer sequences that have been assembled from those reads.

```
cd ~/annotation_course/practical2/RNAseq
mkdir stringtie
cd stringtie

stringtie ../tophat/tophat_out/accepted_hits.bam -o outdir/transcripts.gtf
```

When done you can find your results in the directory ‘outdir’. The file transcripts.gtf includes your assembled transcripts.
As Webapollo doesn't like the gtf format file you should convert it in gff3 format.
```
gxf_to_gff3.pl --gff transcripts.gtf -o transcripts.gff3
``` 
Then, transfer the gff3 file to your computer and load it into [Webapollo](http://annotation-prod.scilifelab.se:8080/NBIS_course/). How well does it compare with your Augustus results? Looking at your results, are you happy with the default values of Stringtie (which we used in this exercise) or is there something you would like to change?

## 2. De-novo transcriptome assembly:

### Trinity

Trinity assemblies can be used as complementary evidence, particularly when trying to polish a gene build with Pasa. Before you start, check how big the raw read data is that you wish to assemble to avoid unreasonably long run times.

```
cd ~/annotation_course/practical2/RNAseq
mkdir trinity
cd trinity

Trinity --seqType fq --max_memory 64G --left ~/annotation_course/data/RNAseq/ERR305399.left.fastq.gz --right ~/annotation_course/data/RNAseq/ERR305399.right.fastq.gz --CPU 8 --output trinity_result --SS_lib_type RF 
```

Trinity takes a long time to run if you want to have a look at the results, look in ~/annotation_course/course_material/data/dmel/chromosome_4/RNAseq/ the output that will be used later on for the annotation will be Trinity.fasta

## Closing remarks

You have now successfully perform transcript assemblies. You have seen how to perform a genome-guided assembly as well as de-no assembly.
