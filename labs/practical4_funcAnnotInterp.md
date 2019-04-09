# Functional annotation

Functional annotation is the process during which we try to put names to faces - what do genes that we have annotated and curated? Basically all existing approaches accomplish this by means of similarity. If a translation product has strong similarity to a protein that has previously been assigned a function, the function in this newly annotated transcript is probably the same. Of course, this thinking is a bit problematic (where do other functional annotations come from...?) and the method will break down the more distant a newly annotated genome is to existing reference data. A complementary strategy is to scan for more limited similarity - specifically to look for the motifs of functionally characterized protein domains. It doesn't directly tell you what the protein is doing exactly, but it can provide some first indication.

In this exercise we will use an approach that combines the search for full-sequence simliarity by means of 'Blast' against large public databases with more targeted characterization of functional elements through the InterproScan pipeline. Interproscan is a meta-search engine that can compare protein queries against numerous databases. The output from Blast and Interproscan can then be used to add some information to our annotation.

## Prepare the input data

Since we do not wish to spend too much time on this, we will again limit our analysis to chromosome 4. It is also probably best to choose the analysis with ab-initio predictions enabled (unless you found the other build to be more convincing). Maker produces a protein fasta file (called "annotations.proteins.fa") together with the annotation and this file should be located in your maker directory.

Move in the proper folder:  
```
mkdir -p ~/annotation_course/practical4
cd ~/annotation_course/practical4
```
Now link the annotation you choose to work with. The command will looks like:
```
ln -s ~/annotation_course/practical2/complement/maker_abinitio_cplt_by_evidence.gff maker_final.gff  
ln -s ~/annotation_course/practical2/complement/maker_abinitio_cplt_by_evidence.fasta maker_final.faa
```
## Interproscan approach
 Interproscan combines a number of searches for conserved motifs and curated data sets of protein clusters etc. This step may take fairly long time. It is recommended to paralellize it for huge amount of data by doing analysis of chunks of tens or hundreds proteins.

### Perform [InterproScan](https://github.com/ebi-pf-team/interproscan/wiki) analysis
InterproScan can be run through a website or from the command line on a linux server. Here we are interested in the command line approach.
<u>Interproscan allows to look up pathways, families, domains, sites, repeats, structural domains and other sequence features.</u>  

Launch Interproscan with the option -h if you want have a look about all the parameters.

- The '-app' option allows defining the database used. Here we will use the PfamA,ProDom and SuperFamily databases.  
- Interproscan uses an internal database that related entries in public databases to established GO terms. By running the '-goterms' option, we can add this information to our data set.
- If you enable the InterPro lookup ('-iprlookup'), you can also get the InterPro identifier corresponding to each motif retrieved: for example, the same motif is known as PF01623 in Pfam and as IPR002568 in InterPro.
- The option '-pa' provides mappings from matches to pathway information (MetaCyc,UniPathway,KEGG,Reactome).
```
interproscan.sh -i maker_final.fa -t p -dp -pa -appl Pfam,ProDom-2006.1,SuperFamily-1.75 --goterms --iprlookup
```
The analysis shoud take 2-3 secs per protein request - depending on how many sequences you have submitted, you can make a fairly deducted guess regarding the running time.  
You will obtain 3 result files with the following extension '.gff3', '.tsv' and '.xml'. Explanation of these output are available [>>here<<](https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats).


### load the retrieved functional information in your annotation file:
Next, you could write scripts of your own to merge interproscan output into your annotation. Incidentially, Maker comes with utility scripts that can take InterProscan output and add it to a Maker annotation file (you need to load maker).  

- ipr\_update\_gff: adds searchable tags to the gene and mRNA features in the GFF3 files.  
- iprscan2gff3: adds physical viewable features for domains that can be displayed in JBrowse, Gbrowse, and Web Apollo.
```
gff3_sp_manage_functional_annotation.pl --gff maker_final.gff -i maker_final.faa.tsv -o  maker_final.interpro
```
Where a match is found, the new file will now include features called Dbxref and/or Ontology_term in the gene and transcript feature field (9th column).
The improved annotation is the gff file inside the maker_final.interpro folder.

## BLAST approach
Blast searches provide an indication about potential homology to known proteins.
A 'full' Blast analysis can run for several days and consume several GB of Ram. Consequently, for a huge amount of data it is recommended to parallelize this step doing analysis of chunks of tens or hundreds proteins. This approach can be used to give a name to the genes and a function to the transcripts.

### Perform Blast searches from the command line on Uppmax:

To run Blast on your data, use the Ncbi Blast+ package against a Drosophila-specific database (included in the folder we have provided for you, under **annotation\_course/data/blastdb/uniprot\_dmel/uniprot\_dmel.fa**) - of course, any other NCBI database would also work:
```
blastp -db ~/annotation_course/data/blastdb/uniprot_dmel/uniprot_dmel.fa -query maker_final.faa -outfmt 6 -out blast.out -num_threads 8
```
Against the Drosophila-specific database, the blast search takes about 2 secs per protein request - depending on how many sequences you have submitted, you can make a fairly deducted guess regarding the running time.

### load the retrieved information in your annotation file:  

Now you should be able to use the following script:
```
gff3_sp_manage_functional_annotation.pl -f maker_final.interpro.gff -b blast.out --db  ~/annotation_course/data/blastdb/uniprot_dmel/uniprot_dmel.faa -o maker_final.interpro.blast  
```
That will add the name attribute to the "gene" feature and the description attribute (corresponding to the product information) to the "mRNA" feature into you annotation file. 
The improved annotation is the gff file inside the maker_final.interpro.blast folder.

 * How many genes do not have any names ?
 
### Set nice IDs

The purpose is to modify the ID value by something more convenient (i.e FLYG00000001 instead of maker-4-exonerate_protein2genome-gene-8.41).  
```
gff3_sp_manage_functional_annotation.pl -f maker_final.interpro.blast/maker_final.gff --ID FLY -o maker_final.interpro.blast.ID  
```
The improved annotation is the gff file inside the maker_final.interpro.blast.ID folder.

### Polish your file for a nice display within Webapollo

For a nice display of a gff file within Webapollo some modification might be needed.
As example the attribute ***product*** is not displayed in Webapollo, whereas renaming it ***description*** will work out.
```
~/annotation_course/GAAS/annotation/WebApollo/gff3_webApollo_compliant.pl -gff maker_final.interpro.blast.ID/maker_final.gff -o final_annotation.gff
```

## Visualise the final annotation

Transfer the final_annotation.gff file to your computer using scp in a new terminal:

scp -i ~/.ssh/azure_rsa student@__IP__:/home/student/annotation_course/practical4/final_annotation.gff .

Load the file in into the genome portal called drosophila_melanogaster_chr4 in the Webapollo genome browser available at the address [http://annotation-prod.scilifelab.se:8080/NBIS_course/](http://annotation-prod.scilifelab.se:8080/NBIS_course/). [Here find the WebApollo instruction](UsingWebapollo)

Wondeful ! insn't it ?

## What's next?

Because of Makers' compatibility with GMOD standards, an annotation augmented in one or both of this way can be loaded into e.g. WebApollo and will save annotators a lot of work when e.g. adding meta data to transcript models.



# Submission to public repository (creation of an EMBL file)

Once your are satisfied by the wonderful annotation you have done, it would useful important to submit it to a public repostiroy. Fisrt you will be applaused by the community because you share your nice work, secondly this is often mandatory if you wish to publish some work related to this annotation.

Current state-of-the-art genome annotation tools use the GFF3 format as output, while this format is not accepted as submission format by the International Nucleotide Sequence Database Collaboration (INSDC) databases. Converting the GFF3 format to a format accepted by one of the three INSDC databases is a key step in the achievement of genome annotation projects. However, the flexibility existing in the GFF3 format makes this conversion task difficult to perform.

In order to submit to **NCBI**, the use of a tool like [GAG](https://genomeannotation.github.io/GAG/) will save you lot time.  
In order to submit to **EBI**, the use of a tool like [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) will be your best choice.

**Let's prepare your annotation to submit to ENA (EBI)**

In real life, prior to a submission to ENA, you need to create an account and create a project asking a locus_tag for your annotation. You have also to fill lot of metada information related to the assembly and so on. We will skip those tasks using fake information.
First you need to download and install EMBLmyGFF3:
```
pip install --user git+https://github.com/NBISweden/EMBLmyGFF3.git
EMBLmyGFF3 finalOutputDir/codingGeneFeatures.gff 4.fa -o my_annotation_ready_to_submit.embl
```

You now have a EMBL flat file ready to submit. In theory to finsish the submission, you will have to send this archived file to their ftp server and finish the submission process in the website side too.
But we will not go further. We are done. CONGRATULATION you know most of the secrets needed to understand the annotations on and perform your own !
