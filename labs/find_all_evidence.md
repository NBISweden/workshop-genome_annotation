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

## 1 Obtaining Protein

**Swissprot:**  
Uniprot is an excellent source for high quality protein sequences. The main site can be found at [http://www.uniprot.org](http://www.uniprot.org). This is also the place to find Swissprot, a collection of manually curated non-redundant proteins that cover a wide range of organisms while still being manageable in size.

**_Exercise 1_ - Swissprot:**  
Navigate the Uniprot site to find the download location for Swissprot in fasta-format. You do not need to download the file, just find it. In what way does Swissprot differ from Uniref (another excellent source of proteins, also available at the same site)?

**Uniprot:**  
Even with Swissprot available, you also often want to include protein sequences from organisms closely related to your study organism. An approach we often use is to concatenate Swissprot with a few protein fasta-files from closely related organisms and use this in our annotation pipeline.

**_Exercise 2_ - Uniprot:**  
Use Uniprot to find (not download) all protein sequences for all the complete genomes in the family Drosophilidae. How many complete genomes in Drosophilidae do you find?

**Refseq:**  
Refseq is another good place to find non-redundant protein sequences to use in your project. The sequences are to some extent sorted by organismal group, but only to very large and inclusive groups. The best way to download large datasets from refseq is using their ftp-server at [ftp://ftp.ncbi.nlm.nih.gov/refseq/](ftp://ftp.ncbi.nlm.nih.gov/refseq/).

**_Exercise 3_ - Refseq:**  
Navigate the Refseq ftp site to find the invertebrate collection of protein sequences. You do not need to download the sequences, just find them. The files are mixed with other types of data, which files include the protein sequences?

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

**_Exercise 6_ - EBI:**
Go to [ENA website](https://www.ebi.ac.uk/ena) and find how many paired illumina HiSeq 4000 runs are available.
<details>
<summary>:key: Click to see how to get the answer.</summary>
click on the &lt;<strong>search and browse</strong>&gt; tab, then under the &lt;<strong>Free text search</strong>&gt; paragraph click on the &lt;<strong>ENA Advanced Search</strong>&gt; link. You should end up on this page: <i>https://www.ebi.ac.uk/ena/data/warehouse/search</i>.  
  From here two solutions:
  <ol>
   <li>Select <strong>Read</strong> from the &lt;select domain&gt; list; Type <strong>Drosophila melanogaster</strong> into the &lt;Drosophila melanogaster&gt; field;  Select <strong>Paired</strong> from the &lt;Librairy layout&gt; field; select <strong>Illumina HiSeq 4000</strong> from the &lt;Instrument model&gt; field; and click on search.
   <li>Write directly into the &lt;Search query&gt; box the following comand:  
  ```
  library_layout="PAIRED" AND tax_eq(7227) AND instrument_model="Illumina HiSeq 4000"
  ```
  </ol>
</details>  
