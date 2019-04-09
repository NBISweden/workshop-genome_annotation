
# Configure the maker_opts.ctl properly for the abinitio evidence-driven annotation:


\#-----Genome (these are always required)  
**genome=genome.fa** #genome sequence (fasta file or fasta embeded in GFF3 file)  
organism\_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

...

\#-----EST Evidence (for best results provide a file for at least one)  
**est=** #set of ESTs or assembled mRNA-seq in fasta format  
altest= #EST/cDNA sequence file in fasta format from an alternate organism  
**est\_gff=stringtie2genome.genome.ok.gff,est2genome.gff** #aligned ESTs or mRNA-seq from an external GFF3 file  
altest\_gff= #aligned ESTs from a closly relate species in GFF3 format

...

\#-----Protein Homology Evidence (for best results provide a file for at least one)  
**protein=** #protein sequence file in fasta format (i.e. from mutiple oransisms)  
**protein\_gff=protein2genome.gff** #aligned protein homology evidence from an external GFF3 file

...

\#-----Repeat Masking (leave values blank to skip repeat masking)  
**model\_org=** #select a model organism for RepBase masking in RepeatMasker  
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker   
**repeat\_protein=** #provide a fasta file of transposable element proteins for RepeatRunner  
**rm\_gff=repeatmasker.genome.gff,repeatrunner.genome.gff** #pre-identified repeat elements from an external GFF3 file  
prok\_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no  
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

...

\#-----Gene Prediction  
snaphmm= #SNAP HMM file  
gmhmm= #GeneMark HMM file  
**augustus\_species=fly** #Augustus gene prediction species model  
fgenesh\_par\_file= #FGENESH parameter file  
pred\_gff= #ab-initio predictions from an external GFF3 file  
model\_gff= #annotated gene models from an external GFF3 file (annotation pass-through)  
**est2genome=0** #infer gene predictions directly from ESTs, 1 = yes, 0 = no  
**protein2genome=0** #infer predictions from protein homology, 1 = yes, 0 = no  
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no  
snoscan\_rrna= #rRNA file to have Snoscan find snoRNAs  
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

...
**keep_preds=1**
...

To better understand the different parameters you can have a look [here](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained) 
