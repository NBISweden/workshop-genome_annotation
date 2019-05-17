---
layout: default-overview
title: Prokaryote annotation
exercises: 60
questions:
  - How to run bacterial annotation?
  - What should I look in my assembly to go forward
objectives:
  - run prokka
---
# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export bacterial_annotation_path=/proj/g2019006/nobackup/$USER/bacterial_annotation
mkdir -p $bacterial_annotation_path
```


# Organizing data

Before going into the exercises below, you should create in your home folder a specific folder for this practical session and copy a folder with the course data using:  

```
cd $bacterial_annotation_path
```

# Bacterial annotation using Prokka

Before running Prokka on genomes assemblies, it is a good step to start with checking the gene content of the assembly

## Checking genes in the assembly

[BUSCO](https://busco.ezlab.org/) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness. Genes that make up the BUSCO sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species.

You will run BUSCO on 3 bacterial assemblies provided (one E coli, one chlamydia and one streptococcus). We will select the lineage set of bacteria.

BUSCO is using augustus to run, as we have no administator rights on uppmax we need to copy the config file of augustus in a folder we can write in and set up the environment.

```
module load BUSCO/3.0.2b

source $BUSCO_SETUP

run_BUSCO.py -i $data/raw_computes/Chlamydia_trachomatis_a_363.fa -o chlamydia_busco -m geno -c 8 -l /sw/apps/bioinfo/BUSCO/v2_lineage_sets/bacteria_odb9
```
look at the results of busco in short_summary_chlamydia_busco.txt

:question:what outputs do you have? where do you have the annotation?
<br>:question:what do you think about this assembly? Is it a good one? can you see any potential problem with continuing the annotation?
<br>:question:how do you expect the annotation will be?

Do the same for the two other assemblies and answer those questions again.

## Prokka

Prokka is a really easy tool to use for bacterial annotation.

You are going to use the same assemblies you used previously for Busco

```
module unload BUSCO

module load prokka/1.12-12547ca

prokka --help
```
The goal of the exercise is for you to learn how to use prokka and to annotate the 3 assemblies and then visualize them in IGV.

run prokka without any options and then with options of your choices (we encourage you to try at least the options --proteins and --rfam )

<details>
<summary>:key: Click to see part of the solution .</summary>  

Running prokka with only the output option looks like this :
<code> prokka $data/raw_computes/Chlamydia_trachomatis_a_363.fa --outdir prokka_Chlamydia
</code>

Running prokka with --proteins and --rfam looks like this :

<code> prokka $data/raw_computes/Chlamydia_trachomatis_a_363.fa --proteins $data/raw_computes/uniprot-chlamydia.fasta --rfam --outdir prokka_Chlamydia_prot_rfam
</code>

You can try other options to see what you would need to modify in your own projects!

</details>


Look at the different results obtained :

:question:Do you see any differences with the different options and no options you used? (like for instance with or without --proteins)
<br>:question:Did you get the annotation you expected after the busco results?

You could now also visualise all this information using a genome browser, such as [IGV](http://software.broadinstitute.org/software/igv/).
IGV requires a genome fasta file and any number of annotation files in GTF or GFF3 format (note that GFF3 formatted file tend to look a bit weird in IGV sometimes).

Transfer the gff3 files to your computer using scp:    
```
scp __YOURLOGIN__@rackham.uppmax.uu.se:/proj/g2019006/nobackup/__YOURLOGIN__/bacterial_annotation/YOURFILE .
```

Congratulations you have annotate bacterial genome!

## Checking gene set completeness (Optional)

BUSCO can also be used after the annotation to check if you found the genes you were expected or if something happened during the annotation and you lost genes. To do so you change the option "-m geno" by "-m prot"

```
module load BUSCO/3.0.2b

run_BUSCO.py -i prokka_Chlamydia/PROKKA_05102019.faa -o chlamydia_busco_prot -m prot -c 8 -l /sw/apps/bioinfo/BUSCO/v2_lineage_sets/bacteria_odb9
```
You can do it for the two other genomes.

:question:Do you see a difference with the BUSCO of the genome?

<details>
<summary>:key: Click to see the solution .</summary>  
Often the BUSCO results for genes are slightly lower than the BUSCO results for the full genome, this is due to the fact that annotation method will always not predict everything.
It should not be too much of a difference either.
Sometimes but really rarely, the BUSCO after annotation will be better than the BUSCO assembly. It is due to the fact that BUSCO check the compleness in two different way. 

</details>
