## Inspect the output

### Finding your way around

By default, Maker will write the output of its different analyses into a folder named:

**&lt;name\_of\_genome\_fasta&gt;.maker.output**

In our case:

**genome.maker.output**

Within the main output directory, Maker keeps a copy of the config files, a database (here: **genome.db**), directories for the blast databases created from your evidence data and a file called **genome\_master\_datastore\_index.log**.

Out of these files, only the **genome\_master\_datastore\_index.log** is really interesting to us. It includes a log of all the contigs included in the genome fasta file - together with their processing status (ideally: FINISHED) and the location of the output files. Since Maker can technically run in parallel on a large number of contigs, it creates separate folders for each of these input data. For larger genomes, this can generate a very deep and confusing folder tree. The **genome\_master\_datastore\_index.log** helps you make sense of it:
```
4       genome_datastore/A8/7F/4/ STARTED  
4       genome_datastore/A8/7F/4/ FINISHED
```
This meens the sequence **4** was started - and finished, with all data (annotation, protein predictions etc) written to the subfolder **genome\_datastore/A8/7F/4/**.

If you look into that folder, you will find the finished Maker annotation for this contig.
```
rw-rw-r- 1 student student 472193 Mar 24 10:16 4.gff
*rw-rw-r- 1 student student 3599 Mar 24 10:16 4.maker.augustus_masked.proteins.fasta
*rw-rw-r- 1 student student 10388 Mar 24 10:16 4.maker.augustus_masked.transcripts.fasta
*rw-rw-r- 1 student student 176 Mar 24 10:16 4.maker.non_overlapping_ab_initio.proteins.fasta 
*rw-rw-r- 1 student student 328 Mar 24 10:16 4.maker.non_overlapping_ab_initio.transcripts.fasta
rw-rw-r- 1 student student 3931 Mar 24 10:16 4.maker.proteins.fasta
rw-rw-r- 1 student student 20865 Mar 24 10:16 4.maker.transcripts.fasta
rw-rw-r- 1 student student 4248 Mar 24 10:15 run.log
drwxrwsr-x 3 student student 4096 Mar 24 10:16 theVoid.4
```
\* only if an abinitio tool has been activated

The main annotation file is '4.gff' - including both the finished gene models and all the raw compute data. The other files include fasta files for the different sequence features that have been annotated - based on ab-initio predictions through augustus as well as on the finished gene models. The folder 'theVoid' include all the raw computations that Maker has performed to synthesize the evidence into gene models.

## Understanding a Maker annotation

You have two options now for gathering the output in some usable form - copy select files by hand to wherever you want them. Or you can use a script that does the job for you (we have included an example in the script folder).

From the folder you have run Maker, run the script called 'maker\_merge\_outputs\_from\_datastore' to create an output file for all annotations and protein files:
```
maker_merge_outputs_from_datastore.pl 
```
This will create a directory called "**maker_output_processed**" containing a long list of files depending on parameters used for running MAKER.  

```
augustus_masked.gff
blastn.gff
blastx.gff
cdna2genome.gff
est_gff_est2genome.gff
evm.gff
fgenesh_masked.gff
genemark.gff
genome.all.maker.augustus_masked.proteins.fasta
genome.all.maker.augustus_masked.transcripts.fasta
genome.all.maker.evm.proteins.fasta
genome.all.maker.evm.transcripts.fasta
genome.all.maker.fgenesh.proteins.fasta
genome.all.maker.fgenesh.transcripts.fasta
genome.all.maker.genemark.proteins.fasta
genome.all.maker.genemark.transcripts.fasta
genome.all.maker.noncoding.fasta
genome.all.maker.non_overlapping_ab_initio.proteins.fasta
genome.all.maker.non_overlapping_ab_initio.transcripts.fasta
genome.all.maker.proteins.fasta
genome.all.maker.snap_masked.proteins.fasta
genome.all.maker.snap_masked.transcripts.fasta
genome.all.maker.transcripts.fasta
genome.all.maker.trnascan.noncoding.fasta
maker_bopts.ctl
maker_evm.ctl
maker_exe.ctl
maker.gff
maker_opts.ctl
maker_stat.txt
protein2genome.gff
protein_gff_protein2genome.gff
repeat_gff_repeatmasker.gff
repeat_gff_repeatrunner.gff
snap_masked.gff
tblastx.gff
```

Here is a describtion of the most important files:

 * **maker.gff** 

It contains the annotation done by maker in ([GFF3 format](http://www.sequenceontology.org/gff3.shtml)). If you use 'less' to read this annotation file, you will see a range of different features:
```
##gff-version 3  
4       maker   gene    24134   25665   .       +       .       ID=maker-4-exonerate_protein2genome-gene-0.0;Name=maker-4-exonerate_protein2genome-gene-0.0
4       maker   mRNA    24134   25665   917     +       .       ID=maker-4-exonerate_protein2genome-gene-0.0-mRNA-1;Parent=maker-4-exonerate_protein2genome-gene-0.0;Name=maker-4-exonerate_protein2genome-gene-0.0-mRNA-1;_AED=0.09;_eAED=0.09;_QI=0|0.33|0.25|1|0|0|4|44|290
```
...

For example, the above lines read:

On the sequence with id ´4´, there is a gene feature located from position 24134 to 25665, on the plus strand and with the id 'maker-4-exonerate\_protein2genome-gene-0.0'. 
On this same sequence, belonging to the gene, is located a transcript from position 24134 to 25665, on the plus strand and with the id 'maker-4-exonerate\_protein2genome-gene-0.0-mRNA-1'. It's quality, or AED score, is 0.09 - which means that the evidence alignments are close to be in perfect agreement with the transcript model.
And so on.

 * **maker_stat.txt**  
Statistics of the annotation (maker.gff annotation file).

 * **genome.all.maker.transcripts.fasta**  
This fasta file contains the nucleotide sequences of the transcripts (mRNA) of the MAKER gene models.

 * **genome.all.maker.proteins.fasta** 
This fasta file contains the amino acid sequences of the proteins translated from the CDS of the MAKER gene models.

 * **cdna2genome.gff**  
Contains the fasta sequences provided by the **est=** parameter that have been succefuly mapped onto the assembly.

 * **est_gff_xxx_.gff**  
Contains the gff data provided by the **est_gff=** parameter that have been kept.

 * **protein2genome.gff**  
Contains the fasta sequences provided by the **protein=** parameter that have been succefuly mapped onto the assembly.

 * **protein_gff_xxx_.gff**  
Contains the gff data provided by the **protein_gff=** parameter that have been kept.

 * **blastx.gff**  
Contains the blast results of the protein fasta sequences provided by the **protein=** parameter.

 * **blastn.gff**  
Contains the blast results of the nucleotide sequences provided by the **est=** parameter.

 * **tblastx.gff**  
Contains the blast results of the nucleotide sequences provided by the **alt_est=** parameter.

 * **repeatmasker.gff**  
Contains the repeat masking results when **model_org=**  or/and  **rmlib=**  parameter is used.

 * **repeatrunner.gff**  
Contains the repeatrunner results when **repeat_protein=** parameter is used.

 * **genome.all.maker.trnascan.noncoding.fasta**  
Contains the result in fasta format of the tRNAscan analysis when the parameter **trna=1** parameter is used. The corresponding gff result is within the maker.gff file.

 * **augustus_masked.gff**  
Contains the raw augustus annotation (not filter by MAKER) un gff format.

 * **genome.all.maker.augustus_masked.transcripts.fasta**  
Contains the transcripts of the raw augustus annotation (not filter by MAKER) in fasta format.

 * **genome.all.maker.augustus_masked.proteins.fasta**  
Contains the proteins of the raw augustus annotation (not filter by MAKER) in fasta format (translated CDS).

 * **maker_xxx.ctl**  
All the file with **.ctl** extensions are copy of the control files used to produce the result contained in the current folder.

