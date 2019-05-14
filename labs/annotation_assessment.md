---
layout: default-overview
title: Comparing and evaluating annotations
exercises: 1h15
questions:
  - How can we compare two annotation?
  - How to merge different annotations?
objectives:
  - learn how to assess the quality of an annotation
  - complement a first annotation with a second one
---

## Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export structural_annotation_path=/proj/g2019006/nobackup/$USER/structural_annotation
export abinitio_augustus_path=/proj/g2019006/nobackup/$USER/abinitio_augustus
```

# Introduction

In this exercise you will handle different annotation files:

• the pure abinitio one done with augustus (practical : Abinitio with augustus),

• the evidence-based done with MAKER

• the abinitio evidence-driven one done with MAKER.

• the official annotation from Ensembl

## overview

Evaluating an annotation can be done in different ways:

• Looking at the number of genes  
It isn't so much a quality check as a measure of congruency - i.e. the resulting numbers don't tell you which of the two gene builds is more correct.

• Comparison with another annotation  
 It does not help neither to see the quality of your annotation but could help to understand the major differences between several annotations.

• Comparison against a reference
 This case is really rare in real life.

• Running busco on proteins obtained from the annotation  
 It provides a nice feeling about the quality of the annotation but is biased by the fact it focus only on well conserved genes between species during evolution. So, what about species specific genes ?

• In reference to the evidence alignments (AED score)  
 It is what Maker uses internally to select gene models. After synthesizing and annotating loci, the resulting model will be ranked against the filtered evidence alignments. The more congruent these two points of information are, the lower the 'annotation edit distance' (AED) will be. The AED score can be used to e.g. check an annotation for problematic models that may then be subjected to manual curation.

### Gene number

As already seen previously you can have a look at the statistics of an annotation with the **gff3_sp_statistics.pl** script.  
As you will note, there are some differences - and of course, this is expected, since different approaches has been used to generate them. The EnsEMBL annotation is originally imported from FlyBase. Obviously, a lot of manual labor and much more data has been put into the FlyBase annotation - and this highlights a common limitation of any computational pipeline. You will simply never reach the same level of quality and detail as seen in a manually curated reference annotation.

### Comparison with another annotation

We will compare the two annotation made with MAKER: the evidence one and the abinitio one.
```
cd $structural_annotation_path/maker
mkdir complement
cd complement
ln -s ../maker_evidence/maker.gff maker_evidence.gff
ln -s ../maker_abinitio/maker.gff maker_abinitio.gff
maker_checkFusionSplitBetweenTwoBuilds.pl --ref maker_evidence.gff --tar maker_abinitio.gff --out maker_evidence_compare_to_abinitio
cat maker_evidence_compare_to_abinitio/resume.txt
```

:question:How many genes are specific to each annotation ?  
:question:How many genes from the evidence annotation have been merged/fused together by the abinitio annotation ?  

Those two annotations have genes that are not in common (non-overlaping). Let's create a non-redundant concatenated gene set:
```
gff3_sp_complement_annotations.pl --ref maker_abinitio.gff --add maker_evidence.gff -o maker_abinitio_cplt_by_evidence.gff
```
:question:How many genes have been added in this new maker_abinitio_cplt_by_evidence.gff annotation ?

Let's extract the proteins form this new annotation:
```
ln -s $data/genome/genome.fa
gff3_sp_extract_sequences.pl -gff maker_abinitio_cplt_by_evidence.gff -f genome.fa -p -o maker_abinitio_cplt_by_evidence.fa
```

### BUSCO

BUSCO is run before annotating to check if the assembly is good and therefore if the annotation will be good. It is also run after the structural annotation to then compare if we indeed find a number of genes corresponding of the first run of busco.

You will need to link the protein file created by maker on the run with the ab-initio
```
cd $structural_annotation_path/maker
mkdir busco
cd busco

ln -s  $structural_annotation_path/maker/complement/maker_abinitio_cplt_by_evidence.fa

module load BUSCO/3.0.2b
source $BUSCO_SETUP

run_BUSCO.py -i maker_abinitio_cplt_by_evidence.fa -o dmel_maker_abinitio_cplt_by_evidence -m prot -c 8 -l /sw/apps/bioinfo/BUSCO/v2_lineage_sets/arthropoda_odb9
```
 :question:if you compare with you first busco results what do you see?

### Comparison with the reference annotation

As with many tasks within bioinformatics, it is always a great idea to first look around for existing solutions. In the case of comparing annotations, there are in fact options already out there. One such example is genometools, which we have briefly used before.  

First create the working folder:
```
cd $structural_annotation_path/maker
mkdir compare_ref
cd compare_ref
```

Then, copy or sym-link the EnsEMBL reference annotation as well as yours:
```
ln -s $abinitio_augustus_pathaugustus_drosophila.gff
ln -s $structural_annotation_path/maker/complement/maker_abinitio_cplt_by_evidence.gff
ln -s $data/annotation/ensembl.genome.gff
```

Now we have to sort any GFF3-formatted annotation in a way that genometools accepts:
```
module load GenomeTools/1.5.9
gt gff3 -sort augustus_drosophila.gff > augustus_drosophila.sorted.gff
gt gff3 -sort maker_abinitio_cplt_by_evidence.gff > maker_abinitio_cplt_by_evidence.sorted.gff
gt gff3 -sort ensembl.genome.gff > ensembl.sorted.gff
```

With the sorted files, we can now perform a comparison two by two:
```
gt eval ensembl.sorted.gff augustus_drosophila.sorted.gff
gt eval ensembl.sorted.gff maker_abinitio_cplt_by_evidence.sorted.gff
```

This will create a long list of measures for all relevant sequence features with respect to both the 'sensitivity' and 'specificity' - as a measure of how close the annotation comes to a reference. As a reminder, 'specificity' measures the fraction of a reference overlapping a prediction whereas 'sensitivity' measures the fraction of a prediction overlapping a reference.

Note that the measures employed by genometools function in a all-or-nothing fashion. If the overlap is not 100%, it doesn't count (which is why you are unlikely to find gene-level congruencies between your gene builds and the reference annotation).  

:question:From the comparison of your annotations (the pure abinitio Augustus one and the one made with MAKER) to the Ensembl annotation, which one **seems** to be the most comprehensive to you ?

### Filter MAKER annotation by AED score

A AED value of 0 means the whole gene model is supported by evidence while 1 means there is none. Let's try to select only models with good congruency with evidence lines, AED <0.3.

```
cd $structural_annotation_path/maker/
mkdir filter
cd filter
ln -s $structural_annotation_path/maker/complement/maker_abinitio_cplt_by_evidence.gff
maker_select_models_by_AED_score.pl -f maker_abinitio_cplt_by_evidence.gff -v 0.3 -t "<" -o result
```

:question:How many genes have passed your filter ? How many have been discarded ?

## Visualising annotations

**Note:** The following section overlaps with some of the exercises you have done earlier (comparing augustus predictions against the reference annotation).

In the previous tasks, we have looked at the overlap between different gene builds. While this gives us an indication of how similar two annotations are, **it doesn't really allow us to judge the overall quality and similarity of annotations**. Remember, sensitivity and specificity are 'all-or-nothing' - two annotations may be considered very different, but provide similar information, biologically. By that, we mean that two gene models don't need to be 100% identical in their coordinates to tell the scientist that a gene indeed exists in a given location and what it's product looks like.

We therefore need to visually inspect and compare the gene builds. This is a crucial step in any annotation project - gene build pipelines use a set of defined rules, but human pattern recognition is needed to spot potential systematic errors. For example, a pipeline like Maker will simply take all your input and try to synthesize it into an annotation, but it doesn't do too much checks on the data itself. What if you RNA-seq data is messier than you thought? What if your protein data set includes to many 'predicted' proteins that are in clear conflict with the other data?

There exist a number of 'annotation viewers' - IGV, Argo and Apollo, to name a few. A common choice for annotators is the web-based version of Apollo, WebApollo, mostly for its curation capabilities.

### Using WebApollo to view annotations
Transfer your maker annotation files to your computer using the scp command.  
Then, jump to [WebApollo](http://annotation-prod.scilifelab.se:8080/NBIS_course/) and upload your annotation track into the genome portal called **drosophila\_melanogaster\_chr4**. [Here find the WebApollo instruction](labs/webapollo_usage.md)  
You can now compare your gene builds against this reference. Some questions to ask yourself:

:question:Do my gene builds recover all the genes found in the reference?  
:question:What sort of differences are most common?  
