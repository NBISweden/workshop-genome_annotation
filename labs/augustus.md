# Running an ab initio gene finder

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](LoginInstructions).

```
cd ~/annotation_course/

mkdir abinitio_augustus

cd abinitio_augustus
```

We have made a genome browser called Webapollo available for you on the address [http://annotation-prod.scilifelab.se:8080/NBIS_gp1/](http://annotation-prod.scilifelab.se:8080/NBIS_gp1/)  called drosophila\_melanogaster\_course.
This browser can already has a number of tracks preloaded for you, but you can also load data you have generated yourself using the ‘file” menu and then ‘open’ and ‘local files’. First time you go there you need to log in using the email adress provided to register the course and your last name as password (lower case and if more than one last name separated by _ eg: lastname1_lastname2)(if you already have access to our webapollo please use the password that have been previously provided to you).

<u>**Ab initio gene finders:**</u> These methods have been around for a very long time, and there are many different programs to try. We will in this exercise focus on the gene finder Augustus. These gene finders use likelihoods to find the most likely genes in the genome. They are aware of start and stop codons and splice sites, and will only try to predict genes that follow these rules. The most important factor here is that the gene finder needs to be trained on the organism you are running the program on, otherwise the probabilities for introns, exons, etc. will not be correct. Luckily, these training files are available for Drosophila.

**_Exercise 1_ - Augustus:**

First you need to be sure that you have access to the libraries required to run tools (you need to redo this if you have been logged off).

```
export PERL5LIB=$PERL5LIB:~/annotation_course/GAAS/annotation/
```


Second load the needed modules using:  
```
module load bioinfo-tools  
module load augustus
```
Run Augustus on your genome file using:  
```
augustus --species=fly ~/annotation_course/data/genome.fa > augustus_drosophila.gtf
```

Take a look at the result file using ‘less augustus\_drosophila.gtf’. What kinds of features have been annotated? Does it tell you anything about UTRs?

The gff-format of Augustus is non-standard (looks like gtf) so to view it in a genome browser you need to convert it. You can do this using genometools which is available on Uppmax.

Do this to convert your Augustus-file:

```
module load perl  
module load perl_modules  
module load BioPerl/1.6.924_Perl5.18.4

~/annotation\_course/course\_material/git/GAAS/annotation/Tools/Converter/gxf_to_gff3.pl -g augustus_drosophila.gtf -o augustus_drosophila.gff3 --gff_version 2
```
Transfer the augustus\_drosophila.gff3 to your computer using scp:    
```
scp __login__@rackham.uppmax.uu.se:~/annotation_course/practical2/augustus_drosophila.gff3 .  
```
Load the file in [Webapollo](http://annotation-prod.scilifelab.se:8080/NBIS_gp1/). [Here find the WebApollo instruction](webapollo_usage)
<br/>Load the Ensembl annotation available in  ~/annotation\_course/course\_material/data/dmel/chromosome\_4/annotation
How does the Augustus annotation compare with the Ensembl annotation? Are they identical?

**_Exercise 2 -_ Augustus with yeast models:**  
Run augustus on the same genome file but using settings for yeast instead (change species to Saccharomyces).

Load this result file into Webapollo and compare with your earlier results. Can you based on this draw any conclusions about how a typical yeast gene differs from a typical Drosophila gene?
