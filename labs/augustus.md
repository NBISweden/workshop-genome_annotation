---
layout: default-overview
title: Abinitio with augustus
exercises: 45
questions:
  - How to run augustus?
  - How visualise you results?
objectives:
  - Run augustus
  - Open and load data in webapollo
---


# Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export abinitio_augustus_path=/proj/g2019006/nobackup/$USER/abinitio_augustus
mkdir -p $abinitio_augustus_path
```

# Running an ab initio gene finder

```bash
cd $abinitio_augustus_path
```

<u><strong>Ab initio gene finders:</strong></u> These methods have been around for a very long time, and there are many different programs to try. We will in this exercise focus on the gene finder Augustus. These gene finders use likelihoods to find the most likely genes in the genome. They are aware of start and stop codons and splice sites, and will only try to predict genes that follow these rules. The most important factor here is that the gene finder needs to be trained on the organism you are running the program on, otherwise the probabilities for introns, exons, etc. will not be correct. Luckily, these training files are available for Drosophila.


:mortar_board: **Augustus:**

First you need to be sure that you have access to the libraries required to run tools (you need to redo this if you have been logged off).

Second load the needed modules using:  
``` bash
module load augustus/3.2.3
```
Then you can have a look at the list of species that already have a trained hmm model.  

```bash
augustus --species=help
```

:question:Did you see the appropriate model for Drosophila Melanogaster?

So, let's now launch Augustus on our genome with the `fly` model.

```bash
augustus --species=fly $data/genome/genome.fa --gff3=yes --progress=true > augustus_drosophila.gff
```

if you wish to annotate isoforms too, use the following command:

```bash
augustus --species=fly $data/genome/genome.fa --gff3=yes --progress=true --alternatives-from-sampling=true > augustus_drosophila_isoform.gff
```

Take a look at the result file using ‘less augustus\_drosophila.gff’. What kinds of features have been annotated? Does it tell you anything about UTRs?

The gff-format of Augustus is non-standard (looks like gtf) so to view it in a genome browser you need to convert it. You can do this using the following command line:

```bash
gxf_to_gff3.pl -g augustus_drosophila.gff -o augustus_drosophila.gff3
```
To better understand what contains your gff file you may use a script that will provide you some statistics like this one:
```bash
gff3_sp_statistics.pl --gff augustus_drosophila.gff
```
:question:How many genes have you annotated?


Transfer the augustus\_drosophila.gff3 to your computer using scp:    
```bash
scp __YOURLOGIN__@rackham.uppmax.uu.se:/proj/g2019006/nobackup/__YOURLOGIN__/abinitio_augustus/augustus_drosophila.gff3 .  
```
Load the file in [Webapollo](http://annotation-prod.scilifelab.se:8080/NBIS_course). [Here find the WebApollo instruction](webapollo_usage)
<br/>The official Ensembl annotation is available in the genome browser.
:question: How does the Augustus annotation compare with the Ensembl annotation? Are they identical?

:mortar_board: **Augustus with yeast models:**  
Run augustus on the same genome file but using settings for yeast instead (change species to Saccharomyces).

<details>
<summary>:key: Click to see the solution .</summary>
<code> augustus --species=saccharomyces $data/genome/genome.fa --gff3=on > augustus_saccharomyces.gff

</code>
</details>

Load this result file into Webapollo and compare with your earlier results.
:question: Can you based on this draw any conclusions about how a typical yeast gene differs from a typical Drosophila gene?

# Closing remarks

We have seen how to assess the quality of the assembly and how to launch a quick annotation using an abinitio tool.
We have also seen the importance to use a species specific hmm model into the ab initio tool. Thus, the limitation of this approach is linked to the pre-trained species that are available.
