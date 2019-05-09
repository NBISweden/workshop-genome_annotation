---
layout: default-overview
title: Find evidence
exercises: 40
questions:
  - What data do I need to run a genome annotation?
  - How/where to get them?
objectives:
  - Navigate in the different database/websites that provide data
---

## 1 Obtaining Protein (Uniprot)
The mission of UniProt is to provide the scientific community with a comprehensive, high-quality and freely accessible resource of protein sequence and functional information.  
The site can be found at [http://www.uniprot.org](http://www.uniprot.org).
Statistics about the different DB holded by Uniprot can be found on this page [https://www.uniprot.org/statistics/](https://www.uniprot.org/statistics/)

**UniProtKB/Swiss-prot:**  
Swiss-Prot (created in 1986) is the manually annotated and reviewed section of the UniProt Knowledgebase (UniProtKB). It is a high quality annotated and non-redundant protein sequence database, which brings together experimental results, computed features and scientific conclusions. Since 2002, it is maintained by the UniProt consortium and is accessible via the UniProt website.  

:question: How many proteins Hare there in UniProtKB/Swiss-prot? Navigate the Uniprot site to find the download location for Swissprot in fasta-format. You do not need to download the file, just find it.

**UniProtKB/TrEMBL**  
UniProtKB/TrEMBL contains the translations of all coding sequences (CDS) present in the EMBL/GenBank/DDBJ Nucleotide Sequence Databases (INSDC consortium) and also protein sequences extracted from the literature or submitted to UniProtKB/Swiss-Prot.  

Even with UniProtKB/Swiss-prot available, you also often want to include protein sequences from organisms closely related to your study organism. An approach we often use is to concatenate Swissprot with a few protein fasta-files from closely related organisms and use this in our annotation pipeline.

:question: How many proteins Hare there in UniProtKB/TrEMBL? Find (not download) all unrewied protein sequences that have atleast an evidence of their existence at transcript level and come from species of the Drosophilidae family.

<details>
<summary>:key: Click to see how to get the answer.</summary>
On the &lt;<strong>search bar</strong>&gt;, click the &lt;<strong>advanced</strong>&gt; button.  
  From the new opened tab:
  <ol>
   <li>Select <strong>Unreviewed</strong> from the list.</li>
   <li>Type the <strong>+</strong> icon to add one more criteria then select <strong>Protein Existence [PE]>Evidence at transcript level</strong> from the list.</li>
   <li>Type the <strong>+</strong> icon to add one more criteria then select <strong>Taxonomy [OC]</strong> from the list and type Drosophilidae.</li>
  </ol>
The <strong>search query</strong> corresponding to this task is the following:  
       <code>taxonomy:drosophilidae existence:"Evidence at transcript level [2]" AND reviewed:no</code>
</details>

**UniParc**
UniParc is a comprehensive and non-redundant database that contains most of the publicly available protein sequences in the world. UniParc avoids redundancy by storing each unique sequence only once.
If you wish to use a massive database Uniparc will be your best choice. But remember that it probably contains more redundant and spurious data!

:question: How many proteins are there in UniParc?

**Proteomes**
This section contains protein sets from fully sequenced genomes.

:question: How many reference proteomes are there in Proteomes? How many are from eukaryota? What is the size of Drosophila melanogaste proteome? How many of its proteins are reviewed?

**Ensembl:**  
The European Ensembl project makes data available for a number of genome projects, in particular vertebrate animals, through their excellent webinterface. This is a good place to find annotations for model organisms as well as download protein sequences and other types of data. They also supply the Biomart interface, which is excellent if you want to download data for a specific region, a specific gene, or create easily parsable file with gene names etc.

**_Exercise 4_ - Ensembl Biomart:** 
Go to Biomart at [http://www.ensembl.org/biomart/martview](http://www.ensembl.org/biomart/martview) and use it to download all protein sequences for chromosome 4 in Drosophila melanogaster. Once you have downloaded the file, use some command line magic to figure out how many sequences are included in the file. Please ask the teachers if you are having problems here.

## 2 Obtaining EST

EST data is not commonly generated anymore, but may become useful for some projects where such data is still available. Examples may include older genomes targeted for re-annotation or genomes with available EST data for closely related species.

The NCBI or EBI websites are the most appropriate places to retrieve such kind of data.

**_Exercise 5_ - NCBI:**  
Go to the NCBI website and find how many ESTs are available for the drosophila melanogaster species.

## 3 Obtaining RNA-seq

Commonly, such data are produced within the project you are working on. Otherwise the most appropriate data could be retrieved on the Sequence Read Archive (SRA) website from the NCBI or the European Nucleotide Archive (ENA) from the EBI.

:mortar_board: **_Exercise 6_ - EBI:**  
Go to [ENA website](https://www.ebi.ac.uk/ena).  

:question: How many paired illumina HiSeq 4000 runs are available.

<details>
<summary>:key: Click to see how to get the answer.</summary>
click on the &lt;<strong>search and browse</strong>&gt; tab, then under the &lt;<strong>Free text search</strong>&gt; paragraph click on the &lt;<strong>ENA Advanced Search</strong>&gt; link. You should end up on this page: <i>https://www.ebi.ac.uk/ena/data/warehouse/search</i>.  
  From here:
  <ol>
   <li>Select <strong>Read</strong> from the &lt;<strong>select domain</strong>&gt; list.</li>
   <li>Type <strong>Drosophila melanogaster</strong> into the &lt;<strong>Taxon name</strong>&gt; field.</li>
   <li>Select <strong>Paired</strong> from the &lt;<strong>Librairy layout</strong>&gt; field.</li>
   <li>Select <strong>Illumina HiSeq 4000</strong> from the &lt;<strong>Instrument model</strong>&gt; field.</li>
   <li>Click on search.</li>
  </ol>
The <strong>search query</strong> corresponding to this task is the following:  
       <code>library_layout="PAIRED" AND tax_eq(7227) AND instrument_model="Illumina HiSeq 4000"</code>
</details>
