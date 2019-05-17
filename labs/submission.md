---
layout: default-overview
title: Submission
exercises: 20
questions:
  - How to prepare annotation file for submission?
  - How to submit an annotation to INSDC database?
objectives:
  - be able to prepare an annotation file for submission
---

# Prerequisites
For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export submission_path=/proj/g2019006/nobackup/$USER/submission
export structural_annotation_path=/proj/g2019006/nobackup/$USER/structural_annotation
mkdir -p $submission_path
cd $submission_path
ln -s $data/genome/genome.fa
ln -s $structural_annotation_path/maker/complement/maker_abinitio_cplt_by_evidence.gff maker_final.gff

```

<u>**Setup:**</u> For this exercise you need to be logged in to Uppmax. Follow the [UPPMAX login instructions](uppmax_login).

# Submission to public repository (creation of an EMBL file)

Once your are satisfied by the wonderful annotation you have done, it would useful important to submit it to a public repostiroy. Fisrt you will be applaused by the community because you share your nice work, secondly this is often mandatory if you wish to publish some work related to this annotation.

Current state-of-the-art genome annotation tools use the GFF3 format as output, while this format is not accepted as submission format by the International Nucleotide Sequence Database Collaboration (INSDC) databases. Converting the GFF3 format to a format accepted by one of the three INSDC databases is a key step in the achievement of genome annotation projects. However, the flexibility existing in the GFF3 format makes this conversion task difficult to perform.

In order to submit to **NCBI**, the use of a tool like [GAG](https://genomeannotation.github.io/GAG/) will save you lot time.  
In order to submit to **EBI**, the use of a tool like [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) will be your best choice.

**Let's prepare your annotation to submit to ENA (EBI)**

In real life, prior to a submission to ENA, you need to create an account and create a project asking a locus_tag for your annotation. You have also to fill lot of metada information related to the assembly and so on. We will skip those tasks using fake information.

First you need to filter and add extra information to your file otherwise submission might fail:
```bash
gff3_sp_flag_short_introns.pl --gff maker_final.gff -o maker_final_short_intron_flagged.gff
gff3_sp_fix_features_locations_duplicated.pl --gff -o maker_final_short_intron_flagged_duplicated_location_fixed.gff
```

Then you need to download and install EMBLmyGFF3:
```bash
module load python/2.7.6
pip install --user git+https://github.com/NBISweden/EMBLmyGFF3.git
~/.local/bin/EMBLmyGFF3 maker_final_short_intron_flagged_duplicated_location_fixed.gff genome.fa -o my_annotation_ready_to_submit.embl
```

Before to try to submit your file, you must check that everything is fine with the official embl-api-validator. You can find it at the [ena repository](https://github.com/enasequence/sequencetools). Download the validator and validate your file.
```bash
wget http://central.maven.org/maven2/uk/ac/ebi/ena/sequence/embl-api-validator/1.1.265/embl-api-validator-1.1.265.jar
java -jar embl-api-validator-1.1.265.jar -r my_annotation_ready_to_submit.embl
```

If the file is validated, you now have a EMBL flat file ready to submit. In theory to finsish the submission, you will have to send this archived file to their ftp server and finish the submission process in the website side too.
But we will not go further. We are done. CONGRATULATION you know most of the secrets needed to understand the annotations on and perform your own !
