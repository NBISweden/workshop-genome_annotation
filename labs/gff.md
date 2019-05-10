---
layout: default-overview
title: gff/gtf
exercises: 30
questions:
  - How to recognize the format?
  - How to fix it?
  - What information could I extract from those files?
objectives:
  - Understand gff and gtf formats
---

# Prerequisites

For this exercise you need to be logged in to Uppmax.

Setup the folder structure:

```bash
source ~/git/GAAS/profiles/activate_rackham_env
export data=/proj/g2019006/nobackup/$USER/data
export work_with_gff=/proj/g2019006/nobackup/$USER/work_with_gff
mkdir -p $work_with_gff
```

# Recognize the format

```bash
cd $work_with_gff
cp $data/annotation/augustus.xxx .
nano augustus.xxx
```

:question: What is the format of this file? Do you see any problem?

<details>
<summary>:key: Click to see the solution .</summary>
<code> augustus --species=saccharomyces $data/genome/genome.fa --gff3=on > augustus_saccharomyces.gff

</code>
</details>
