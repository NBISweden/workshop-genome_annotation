# Making an abinitio evidence-driven annotation with MAKER

The recommended way of running Maker is in combination with one or more ab-initio profile models. Maker natively supports input from several tools, including augustus, snap and genemark. The choice of tool depends a bit on the organism that you are annotating - for example, GeneMark-ES is mostly recommended for fungi, whereas augustus and snap have a more general use.

The biggest problem with ab-initio models is the process of training them. It is usually recommended to have somewhere around 500-1000 curated gene models for this purpose. Naturally, this is a bit of a contradiction for a not-yet annotated genome.

However, if one or more good ab-initio profiles are available, they can potentially greatly enhance the quality of an annotation by filling in the blanks left by missing evidence. Interestingly, Maker even works with ab-initio profiles from somewhat distantly related species since it can create so-called hints from the evidence alignments, which the gene predictor can take into account to fine-tune the predictions.

Usually when no close ab-initio profile exists for the investigated species, we use the first round of annotation (evidence based) to create one. We first filter the best gene models from this annotation, which are used then to train the abinitio tools of our choice.

In order to compare the performance of Maker with and without ab-initio predictions in a real-world scenario, we have first run a gene build without ab-initio predictions. Now, we run a similar analysis but enable ab-initio predictions through augustus.

## Prepare the input data

No need to re-compute the mapping/alignment of the different lines of evidence. Indeed, this time consuming task has already been performed during the previous round of annotation (evidence based). So, we will use the corresponding gff files previously produced by MAKER.

Link the gff files you want to use into your folder:

 - repeatmasker.chr4.gff (already present)
 - repeatrunner.chr4.gff (already present)
 - genome.fa (already present)
 - stringtie2genome.genome.gff (already present) 
 - est2genome.gff 
 - protein2genome.gff 

```
ln -s maker_evidence/est2genome.gff 
ln -s maker_evidence/protein2genome.gff
```

This time, we do specify a reference species to be used by augustus, which will enable ab-initio gene finding and keep_preds=1 will also show abinitio prediction not supported by any evidences :  

*augustus\_species=fly* #Augustus gene prediction species model  (this is where you can call the database you trained for augustus)   
...
*keep_preds=1*

We must deactivate the evidence base predidction to enable MAKER to pass those alignalignents/hints to the ab-initio tool (That enable the ab-initio evidence-driven mode. Otherwise it would be pure abinitio).  

<i>protein2genome=0</i>  
<i>est2genome=0</i>


With these settings, Maker will run augustus to predict gene loci, but inform these predictions with information from the protein and est alignments.

Before running MAKER you can check you have modified the maker_opts.ctl file properly [here](practical2_supl3_maker.md).

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
