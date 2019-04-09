# Obtaining Protein

<u>**Swissprot:**</u> Uniprot is an excellent source for high quality protein sequences. The main site can be found at [http://www.uniprot.org](http://www.uniprot.org). This is also the place to find Swissprot, a collection of manually curated non-redundant proteins that cover a wide range of organisms while still being manageable in size.

**_Exercise 1_ - Swissprot:**  
Navigate the Uniprot site to find the download location for Swissprot in fasta-format. You do not need to download the file, just find it. In what way does Swissprot differ from Uniref (another excellent source of proteins, also available at the same site)?

<u>**Uniprot:**</u> Even with Swissprot available, you also often want to include protein sequences from organisms closely related to your study organism. An approach we often use is to concatenate Swissprot with a few protein fasta-files from closely related organisms and use this in our annotation pipeline.

**_Exercise 2_ - Uniprot:**  
Use Uniprot to find (not download) all protein sequences for all the complete genomes in the family Drosophilidae. How many complete genomes in Drosophilidae do you find?

<u>**Refseq:**</u> Refseq is another good place to find non-redundant protein sequences to use in your project. The sequences are to some extent sorted by organismal group, but only to very large and inclusive groups. The best way to download large datasets from refseq is using their ftp-server at [ftp://ftp.ncbi.nlm.nih.gov/refseq/](ftp://ftp.ncbi.nlm.nih.gov/refseq/).

**_Exercise 3_ - Refseq:**  
Navigate the Refseq ftp site to find the invertebrate collection of protein sequences. You do not need to download the sequences, just find them. The files are mixed with other types of data, which files include the protein sequences?

<u>**Ensembl:**</u> The European Ensembl project makes data available for a number of genome projects, in particular vertebrate animals, through their excellent webinterface. This is a good place to find annotations for model organisms as well as download protein sequences and other types of data. They also supply the Biomart interface, which is excellent if you want to download data for a specific region, a specific gene, or create easily parsable file with gene names etc.

**_Exercise 4_ - Ensembl Biomart:**  
Go to Biomart at [http://www.ensembl.org/biomart/martview](http://www.ensembl.org/biomart/martview) and use it to download all protein sequences for chromosome 4 in Drosophila melanogaster. Once you have downloaded the file, use some command line magic to figure out how many sequences are included in the file. Please ask the teachers if you are having problems here.

# Obtaining EST

EST data is not commonly generated anymore, but may become useful for some projects where such data is still available. Examples may include older genomes targeted for re-annotation or genomes with available EST data for closely related species.

The NCBI or EBI websites are the most appropriate places to retrieve such kind of data.

**_Exercise 5_ - NCBI:**  
Go to the NCBI website and find how many ESTs are available for the drosophila melanogaster species.

# Obtaining RNA-seq

Commonly, such data are produced within the project you are working on. Otherwise the most appropriate data could be retrieved on the Sequence Read Archive (SRA) website from the NCBI or the European Nucleotide Archive (ENA) from the EBI.
