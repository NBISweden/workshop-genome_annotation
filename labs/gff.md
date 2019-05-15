---
layout: default-overview
title: gff/gtf
exercises: 30
questions:
  - How to recognize the format?
  - How to fix it?
  - What information could I extract from those files?
objectives:
  - Understand gff and gtf formats
---

# Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export work_with_gff=/proj/g2019006/nobackup/$USER/work_with_gff
mkdir -p $work_with_gff
```

# Recognize the format

GFF and GTF format are close and could be difficult to differentiate. To a complete overview of the format you can have a look in the <strong>cheat sheet</strong> section.

```bash
cd $work_with_gff
cp $data/annotation/augustus.xxx .
less augustus.xxx
```

:question:  
<ol>
   <li>Is it a GFF of GTF file?  </li>
   <li>Do you see any problem in the 3rd colum? </li>
   <li>Which version of the format it is? </li>
   <li>Do you see any problem in the 9th colum? </li>
</ol>

<details>
<summary>:key: Click to see the solution .</summary>
<ol>
<li>This is a <strong>GTF</strong> format. You can see that last column where tag and value are separated by a space (would be a '=' in gf format). Another detail that could help it's the last semi-colon that does not exist within gff format. </li>
<li><strong>gene</strong> and <strong>transcript</strong> are features allowed only in <strong>GTF2.5</strong> while <strong>intron</strong> feature exists only in <strong>GTF1</strong>. <strong>tss</strong> feature do not exist officialy in any version. </li>
<li>Tricky question, it looks like GTF2.5 but it's actually a flavor specific to augustus. </li>
<li>The <strong>gene</strong> and <strong>transcript</strong> features have wrong <strong>attributes</strong>. It is missing the <strong>tag</strong>, they only contain the value. It is suppose to look like <code>tag value</code> </li>
</ol>
</details>  
  
   
Now edit the file to fix the 9th column:  

  ```bash
  chmod +w augustus.xxx
  nano augustus.xxx
  ```

<details>
<summary>:key: Click to see the solution .</summary>
  The two first line must be like that:  
  <code>    
    4	AUGUSTUS        gene    386     13142   0.01    +	.	gene_id g1;<br>
    4	AUGUSTUS        transcript	386     13142   0.01    +	.	transcript_id g1.t1;
  </code>
</details>

Now your file has at least a correct structure!  
Let's convert it to **GFF3** format:  

  ```bash
  gxf_to_gff3.pl --gff augustus.xxx -o augustus.gff3 
  ```

The script **gxf_to_gff3.pl** can be your friend when dealing with GFF/GTF format files. It can deal with any kind of GFF/GTF format (even mixed formats) and errors. It allows to create a standardized **GFF3** format file.

# Extract information from a GFF file

The GFF fomat has been developed to be easy to parse and process by a variety of programs in different languages (e.g Unix tools as grep and sort, perl, awk, etc). For these reasons, they decided that each feature is described on a single line.

Download **human** gff annotation v96 from Ensembl:

```bash
 wget ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz
 ```
 
 :question: What is the size of this file?
 
 Now uncompress it:
 ```bash
 gunzip Homo_sapiens.GRCh38.96.chr.gff3.gz 
 ```
 
 :question: What is the size of the uncompressed file?  
  
 The gff/gtf format has a good compression ratio.
 
 Let's now compute some statistics on this file.
 
 :question:  
<ol>
   <li>How many line are there? </li>
   <li>How many <strong>gene</strong> are there? </li>
   <li>How many <strong>mRNA</strong> are there? </li>
   <li>How many <strong>gene</strong> are there on chrmosome  <strong>1</strong>? </li>
   <li>How many types of feature (3rd column) are there? </li>
</ol>
 
<details>
<summary>:key: Click to see the solution .</summary>
<ol>
<li> <code>wc -l Homo_sapiens.GRCh38.96.chr.gff3</code> </li>
<li> <code>awk '{if($3=="gene") print $0}' Homo_sapiens.GRCh38.96.chr.gff3 | wc -l</code> </li>
<li> <code>awk '{if($3=="mRNA") print $0}' Homo_sapiens.GRCh38.96.chr.gff3 | wc -l</code> </li>
<li> <code>awk '{if($3=="gene" && $1=="1") print $0}' Homo_sapiens.GRCh38.96.chr.gff3 | wc -l</code> </li>
<li> <code>awk '{if($0 !~ /^#/)print $3}' Homo_sapiens.GRCh38.96.chr.gff3 | sort -u</code> </li>
</ol>
</details> 
