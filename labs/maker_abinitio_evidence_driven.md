---
layout: default-overview
title: Making an abinitio evidence-driven annotation with MAKER
exercises: 1h30
questions:
  - How to create structural annotation with evidence and abinitio?
objectives:
  - run maker with augustus
  - understand the parameters files and the output
---

## Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export structural_annotation_path=/proj/g2019006/nobackup/$USER/structural_annotation
cd $structural_annotation_path/maker
module load maker/3.01.2-beta
```

# Introduction

The recommended way of running Maker is in combination with one or more ab-initio profile models. Maker natively supports input from several tools, including augustus, snap and genemark. The choice of tool depends a bit on the organism that you are annotating - for example, GeneMark-ES is mostly recommended for fungi, whereas augustus and snap have a more general use.

The biggest problem with ab-initio models is the process of training them. It is usually recommended to have somewhere around 500-1000 curated gene models for this purpose. Naturally, this is a bit of a contradiction for a not-yet annotated genome.

However, if one or more good ab-initio profiles are available, they can potentially greatly enhance the quality of an annotation by filling in the blanks left by missing evidence. Interestingly, Maker even works with ab-initio profiles from somewhat distantly related species since it can create so-called hints from the evidence alignments, which the gene predictor can take into account to fine-tune the predictions.

Usually when no close ab-initio profile exists for the investigated species, we use the first round of annotation (evidence based) to create one. We first filter the best gene models from this annotation, which are used then to train the abinitio tools of our choice.

In order to compare the performance of Maker with and without ab-initio predictions in a real-world scenario, we have first run a gene build without ab-initio predictions. Now, we run a similar analysis but enable ab-initio predictions through augustus.

## Prepare the input data

No need to re-compute the mapping/alignment of the different lines of evidence. Indeed, this time consuming task has already been performed during the previous round of annotation (evidence based). So, we will use the corresponding gff files previously produced by MAKER.

Link the gff files you want to use into your folder:

 - repeatmasker.genome.gff (already present)
 - repeatrunner.genome.gff (already present)
 - genome.fa (already present)
 - stringtie2genome.gff (already present)
 - est2genome.gff
 - protein2genome.gff

```
ln -s maker_evidence/est2genome.gff
ln -s maker_evidence/protein2genome.gff
```

This time, we do specify a reference species to be used by augustus, which will enable ab-initio gene finding and keep_preds=1 will also show abinitio prediction not supported by any evidences :  

*augustus\_species=fly* #Augustus gene prediction species model  (:bulb:this is where you can call the database you trained for augustus dmel_$USER )

To be able to use dmel_$USER you need to type in the terminal :
```
AUGUSTUS_CONFIG_PATH=/proj/g2019006/nobackup/$USER/structural_annotation/augustus_training/maker_results_noAbinitio_clean/augustus_config
```

...
*keep_preds=1*

We must deactivate the evidence base predidction to enable MAKER to pass those alignaligments/hints to the ab-initio tool (that enables the ab-initio evidence-driven mode. Otherwise it would be pure abinitio).  

<i>protein2genome=0</i>  
<i>est2genome=0</i>


With these settings, Maker will run augustus to predict gene loci, but inform these predictions with information from the protein and est alignments.

Before running MAKER check you have modified the maker_opts.ctl file properly.
<details>
<summary>:key: Click here to see the expected maker_opts.ctl.</summary>
{% highlight bash %}

#-----Genome (these are always required)  
genome=genome.fa #genome sequence (fasta file or fasta embeded in GFF3 file)  
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

...

#-----EST Evidence (for best results provide a file for at least one)  
est= #set of ESTs or assembled mRNA-seq in fasta format  
altest= #EST/cDNA sequence file in fasta format from an alternate organism  
est_gff=stringtie2genome.ok.gff,est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file  
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

...

#-----Protein Homology Evidence (for best results provide a file for at least one)  
protein= #protein sequence file in fasta format (i.e. from mutiple oransisms)  
protein_gff=protein2genome.gff #aligned protein homology evidence from an external GFF3 file

...

#-----Repeat Masking (leave values blank to skip repeat masking)  
model_org= #select a model organism for RepBase masking in RepeatMasker  
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker   
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner  
rm_gff=repeatmasker.genome.gff,repeatrunner.genome.gff** #pre-identified repeat elements from an external GFF3 file  
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no  
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

...

#-----Gene Prediction  
snaphmm= #SNAP HMM file  
gmhmm= #GeneMark HMM file  
augustus_species=fly #Augustus gene prediction species model  
fgenesh_par_file= #FGENESH parameter file  
pred_gff= #ab-initio predictions from an external GFF3 file  
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)  
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no  
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no  
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no  
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs  
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

...
keep_preds=1
...

{% endhighlight %}

</details>  
To better understand the different parameters you can have a look [here](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained)


## Run Maker with ab-initio predictions

With everything configured, run Maker as you did for the previous analysis:
```
maker -c 8
```
We probably expect this to take a little bit longer than before, since we have added another step to our analysis.

Once the run is finished, check that everything went properly. If problems are detected, launch MAKER again.  
```
maker_check_progress.sh
```

## Compile the output

When Maker has finished, compile the output:
```
maker_merge_outputs_from_datastore.pl --output maker_abinitio
```
And again, it is probably best to link the resulting output (maker.gff) to a result folder (the same as defined in the previous exercise e.g. dmel\_results), under a descriptive name.

## Inspect the gene models

To get some statistics of your annotation you could launch :
```
gff3_sp_statistics.pl --gff maker_abinitio/maker.gff
```

We could now also visualise the annotation in the Webapollo genome browser.
