---
layout: default-overview
title: Making an evidence based annotation with MAKER
exercises: 45
questions:
  - how to create a structural annotation based on evidence only?
  -
objectives:
  - Understand maker parameters files
  - run maker
---


## Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export structural_annotation_path=/proj/g2019006/nobackup/$USER/structural_annotation
export RNAseq_assembly_path=/proj/g2019006/nobackup/$USER/RNAseq_assembly
mkdir -p $structural_annotation_path
```

## Overview

The first run of Maker will be done without ab-initio predictions. What are your expectations for the resulting gene build? In essence, we are attempting a purely evidence-based annotation, where the best protein- and EST-alignments are chosen to build the most likely gene models. The purpose of an evidence-based annotation is simple. Basically, you may try to annotate an organism where no usable ab-initio model is available. The evidence-based annotation can then be used to create a set of genes on which a new model could be trained on (using e.g. Snap or Augustus). Selection of genes for training can be based on the annotation edit distance (AED score), which says something about how great the distance between a gene model and the evidence alignments is. A score of 0.0 would essentially say that the final model is in perfect agreement with the evidence.

Let's do this step-by-step:

## Prepare the folder and input data

Create the folder where we will launch this maker run.

```
cd $structural_annotation_path
mkdir maker
cd maker
```

Link the raw computes you want to use into your folder. The files you will need are:

- the gff file of the pre-computed repeats (coordinates of repeatmasked regions)

```
ln -s $data/raw_computes/repeatmasker.genome.gff
ln -s $data/raw_computes/repeatrunner.genome.gff
```

In addition, you will also need the genome sequence.
```
ln -s $data/genome/genome.fa
```
Then you will also need EST and protein fasta file:  
```
ln -s $data/evidence/est.genome.fa
ln -s $data/evidence/proteins.genome.fa
```
To finish you could use a transcriptome assembly (This is the same as the one you made using Stringtie and/or the Trinity.fasta):
```
ln -s $data/RNAseq/stringtie/transcript_stringtie.gff3 stringtie2genome.gff
```
OR
```
ln -s $RNAseq_assembly_path/stringtie/transcript_stringtie.gff3 stringtie2genome.gff
```

/!\\ Always check that the gff files you provides as protein or EST contains match/match_part (gff alignment type ) feature types rather than genes/transcripts (gff annotation type) otherwise MAKER will not use the contained data properly. Here we have to fix the stringtie gff file.

```
gff3_sp_alignment_output_style.pl --gff stringtie2genome.gff -o stringtie2genome.ok.gff
```

You should now have 2 repeat files, 1 EST file, 1 protein file, 1 transcript file, and the genome sequence in the working directory.

For Maker to use this information, we need create the three config files, typing this command:
```
module load maker/3.01.2-beta
maker -CTL
```

You can leave the two files controlling external software behaviors untouched but you need to provide the proper parameters in the file called **maker_opts.ctl**.
To edit the **maker_opts.ctl** file you can use the nano text editor:  
```
nano maker_opts.ctl
```

In the **maker_opts.ctl** you will set:

- name of the genome sequence (genome=)  
- name of the 'EST' file in fasta format  (est=) :bulb:You can write the result of your denovo assembly Trinity.fasta there  
- name of the 'Transcript' file in gff format (est_gff=)  
- name of the 'Protein' set file(s) (protein=)  
- name of the repeatmasker and repeatrunner files (rm_gff=)  

You can list multiple files in one field by separating their names by a **comma** ','.

This time, we do not specify a reference species to be used by augustus, which will disable ab-initio gene finding. Instead we set:

  <i>protein2genome=1</i>  
  <i>est2genome=1</i>

This will enable gene building directly from the evidence alignments.

**/!\ Be sure to have deactivated the parameters _model\_org= #_ and _repeat\_protein= #_ to avoid the heavy work of repeatmasker.**
Before running MAKER check you have modified the maker_opts.ctl file properly.
<details>
<summary>:key: Click here to see the expected maker_opts.ctl.</summary>
{% highlight bash %}

#-----Genome (these are always required)  
genome=genome.fa #genome sequence (fasta file or fasta embeded in GFF3 file)  
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

...

#-----EST Evidence (for best results provide a file for at least one)  
est=est.genome.fa,Trinity.fasta #set of ESTs or assembled mRNA-seq in fasta format  
altest= #EST/cDNA sequence file in fasta format from an alternate organism  
est_gff=stringtie2genome.ok.gff #aligned ESTs or mRNA-seq from an external GFF3 file  
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

...

#-----Protein Homology Evidence (for best results provide a file for at least one)  
protein=proteins.genome.fa #protein sequence file in fasta format (i.e. from mutiple oransisms)  
protein_gff= #aligned protein homology evidence from an external GFF3 file

...

#-----Repeat Masking (leave values blank to skip repeat masking)<br/>
model_org= #select a model organism for RepBase masking in RepeatMasker  
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker   
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner  
rm_gff=repeatmasker.genome.gff,repeatrunner.genome.gff #pre-identified repeat elements from an external GFF3 file  
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no  
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

...

#-----Gene Prediction  
snaphmm= #SNAP HMM file  
gmhmm= #GeneMark HMM file  
augustus_species= #Augustus gene prediction species model  
fgenesh_par_file= #FGENESH parameter file  
pred_gff= #ab-initio predictions from an external GFF3 file  
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)  
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no  
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no  
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no  
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs  
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

...
{% endhighlight %}
</details>  
To better understand the different parameters you can have a look [here](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained)

## Run Maker

If your maker\_opts.ctl is configured correctly, you should be able to run maker:
```
maker -c 8
```
This will start Maker on 8 cores, if everything is configured correctly.
This will take a little while and process a lot of output to the screen. Luckily, much of the heavy work - such as repeat masking - are already done, so the total running time is quite manageable, even on a small number of cores.

Once the run is finished, check that everything went properly. If problems are detected, launch MAKER again.
```
maker_check_progress.sh
```

## Inspect the output (optional)

[Here you can find details about the MAKER output.](maker_output_details.md)

## Compile the output

Once Maker is finished, compile the annotation:
```
maker_merge_outputs_from_datastore.pl --output maker_evidence
```
We have specified a name for the output directory since we will be creating more than one annotation and need to be able to tell them apart.  

This should create a **maker\_evidence** folder containing all computed data including **maker.gff** which is the maker annotation file and **genome.all.maker.proteins.fasta** which is the protein fasta file of this annotation. Those two files are the most important outputs from this analysis.

=> You could sym-link the **maker.gff** and **genome.all.maker.proteins.fasta** files to another folder called e.g. dmel\_results, so everything is in the same place in the end. Just make sure to call the links with specific names, since any maker output will be called similarly.


## Inspect the gene models

To get some statistics of your annotation you could read the **maker_stat.txt** file from the **maker\_evidence** folder or launch this script that work on any gff file :
```
gff3_sp_statistics.pl --gff maker_evidence/maker.gff
```

We could now also visualise the annotation in the Webapollo genome browser.
