# Structural annotation with the MAKER pipeline

## Prerequisites

  * **Connection to your virtual machine**  
Before going into the exercises below you need to connect to your virtual machine Ubuntu 16.04 following the instruction we will provide you.

  * **Create the folder structure**  
Now create and move into the **practical2** folder and you are ready to start !
```
mkdir -p ~/annotation_course/practical2
cd ~/annotation_course/practical2
```

  * **List of tools needed. For your convenience they all hae been pre-installed.**  

    * MAKER
    * augustsus
    * GAAS repository

## Overview

**MAKER** is a computational pipeline to automatically generate annotations from a range of input data - including proteins, ESTs, RNA-seq transcripts and ab-initio gene predictions. During this exercise, you will learn how to use Maker with different forms of input data, and how to judge the quality of the resulting annotations.

The Maker pipeline can work with any combination of the following data sets:

* Proteins from the same species or related species  

* Proteins from more distantly related organisms (e.g. Uniprot/Swissprot)  

* Transcriptome sequences from the same species or very closely related species  

* Ab-initio predictions from one or more tools (directly supported are: Augustus, Snap, GeneMark, Fgenesh)  

At minimum, most annotation projects will run with a protein data set, possibly complemented by some RNA-seq data. Popular examples of this are most of the traditional model systems, including human. However, a potential shortcoming of such approaches is that the comprehensiveness of the annotation depends directly on the input data. This can become a problem if our genome of interest is taxonomically distant to well-sequenced taxonomic groups so that only few protein matches can be found. Likewise, not all genes will be expressed at all times, making the generation of a comprehensive RNA-seq data set for annotation challenging.

We will therefore first run our annotation project in the traditional way, with proteins and ESTs, and then repeat the process with a well-trained ab-initio gene predictor. You can then compare the output to get an idea of how crucial the use of a gene predictor is. However, before we get our hands dirty, we need to understand Maker a little better...

Maker strings together a range of different tools into a complex pipeline (e.g. blast, exonerate, repeatmasker, augustus...), fortunately all its various dependencies have been already installed for you. 
Check that everything is running smoothly by creating the MAKER config files:

```
mkdir -p ~/annotation_course/practical2/maker
cd ~/annotation_course/practical2/maker
maker -CTL
```

## Understanding Makers control files

Makers behaviour and information on input data are specified in one of three control files. These are:

- maker_opts.ctl  
- maker_bopts.ctl  
- maker_exe.ctl

What are these files for?

'maker_exe.ctl' holds information on the location of the various binaries required by Maker (including Blast, Repeatmasker etc). Normally, all information in this file will be extracted from $PATH, so if everything is set up correctly, you will never have to look into this file.

Next, 'maker_bopts.ctl' provides access to a number of settings that control the behaviour of evidence aligners (blast, exonerate). The default settings will usually be fine, but if you want to try to annotate species with greater taxonomic distance to well-sequenced species, it may become necessary to decrease stringency of the e.g. blast alignments.

Finally, 'maker_opts.ctl' holds information on the location of input files and some of the parameters controlling the decision making during the gene building.

## Running Maker - Drosophila genome

We will annotate the genome of the fruit fly _Drosophila melanogaster_. First we will perforn a pure evidence based annotation (without ab-initio predictions) and afterwards with ab-initio.

### 1. Creating an evidence based annotation

[Running Maker with only evidence data](practical2_sub_makerNoAbinit.md)

### 2. Creating an abinition evidence-driven annotation

[Running Maker with ab-initio predictions](practical2_sub_makerAbinit.md)

### 3. Inspecting the output

The running of an annotation pipeline like Maker is not actually very hard. But the complicated work is only beginning. How to we best inspect the gene builds? Count features? Visualize it? Most importantly, what steps do we need to take to create a 'finished' annotation that we can use for scientific analyses?

[Comparing and evaluating annotations](practical2_sub_makerCompareAnnot.md)

## Closing remarks

This concludes the gene building part. We have learned how to use the Maker annotation pipeline and have created gene builds with and without ab-initio predictions. Moreover, we have employed some measures to describe and judge these annotations. An essential part that we decided to leave out is the training of ab-initio gene finders. The reason for this omission was that there isn't really any one best way to do this and your mileage may vary a lot based on your organism and input data. Perhaps the most direct approach available at the moment is a combination of evidence-based annotation with Maker and to use the resulting, crude gene models to train SNAP. Since Maker can improve ab-initio predictions 'on the fly', it can tolerate a bit of noise from a less-than-perfect ab-initio profile. If you are setting out on an annotation project, the NBIS annotation service would be happy to discuss the best approach for your data with you.

With that being said, generating a gene build is only one part of an annotation project. Next, we will inspect the annotation in genome browser and make an attempt at functional inference for the predicted gene models.
