---
layout: default
title:  'Prokka exercise'
---

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](files/LoginInstructions).


# Organizing data

Before going into the exercises below, you should create in your home folder a specific folder for this practical session and copy a folder with the course data using:  

```
mkdir -p ~/annotation_course/practical4
cd ~/annotation_course
```

# Bacterial annotation using Prokka

Before running Prokka on genomes assemblies, it is a good step to start with checking the gene content of the assembly

## Checking genes in the assembly

BUSCO2 provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness. Genes that make up the BUSCO2 sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species. It includes 1,066 genes for arthropods, 2,586 for vertebrates, 978 for metazoans, 290 for fungi, 303 for eukaryotes and for bacteria 40 universal marker genes.

You will run BUSCO on 3 bacterial assemblies provided (one E coli, one chlamydia and one streptococcus). We will select the lineage set of bacteria.

BUSCO2 is using augustus to run, as we have no administator rights on uppmax we need to copy the config file of augustus in a folder we can write in and set up the environment.

```
cp -r ~/annotation_course/course_material/augustus_path .*

chmod ug+w -R augustus_path

module load bioinfo-tools 
module load BUSCO 

AUGUSTUS_CONFIG_PATH=augustus_path*

BUSCO -i /home/__login__/annotation\_course/course\_material/data/prokka/Chlamydia_trachomatis_a_363.fa -o chlamydia_busco -m geno -c 8 -l /sw/apps/bioinfo/BUSCO/v2_lineage_sets/bacteria_odb9
```
look at the results of busco in short_summary_chlamydia_busco.txt

- what do you see? 
- what do you think about this assembly? Is it a good one? can you see any potential problem with continuing the annotation?
- how do you expect the annotation will be?

Do the same for the two other assemblies and answer those questions again.

## Prokka

Prokka is a really easy tool to use for bacterial annotation.

You are going to use the same assemblies you used previously for Busco

```
module load bioinfo-tools  
module unload BUSCO

module load prokka

prokka --help
```
The goal of the exercise is for you to learn how to use prokka and to annotate the 3 assemblies and then visualize them in IGV.

run prokka without any options and then with options of your choices (we encourage you to try at least the options --proteins and --rfam )

Look at the different results obtained :

- Do you see any differences with the different options and no options you used?
- Did you get the annotation you expected after the busco results?

You could now also visualise all this information using a genome browser, such as [IGV](http://software.broadinstitute.org/software/igv/). 
IGV requires a genome fasta file and any number of annotation files in GTF or GFF3 format (note that GFF3 formatted file tend to look a bit weird in IGV sometimes).

Transfer the gff3 files to your computer using scp:    
```
scp __login__@milou.uppmax.uu.se:/home/__login__/annotation\_course/practical1/prokka/YOURFILE .
```
- Do you see any differences with the different options and no option you used?

