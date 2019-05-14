---
layout: default-overview
title: Get started
exercises: 20
questions:
  - How do I connect to UPPMAX?
  - How do I get the data and tools necessary for the course?
objectives:
  - Connect to UPPMAX
  - Get all the data necessary and the tools running for the course
---

## Foreword:

We will for all exercises use data for the fruit fly, *Drosophila melanogaster*, as that is one of the currently best annotated organisms and there is plenty of high quality data available. However, working on eukaryotes can be time consuming. Even a small genome like Drosophila would take too long to run within the time we have for this course. Thus to be sure to perform the practicals in good conditions, we will use the smallest chromosome of the drosophila (chromosome 4) like it was a whole genome.

## Prerequisites

### Connection to Uppmax  
Please connect yourself to Uppmax following those instruction [UPPMAX login instructions](../uppmax_login).

### Install dependencies (the GAAS repository through dedicated environment)  
Once connected you will need you to clone and install our github repository indeed this repository contains NBIS annotation team scripts and will be used during this pratical.

  ```bash
  mkdir -p ~/git ; cd ~/git
  git clone https://github.com/NBISweden/GAAS.git
  cd GAAS
  make install
  source ~/git/GAAS/profiles/activate_rackham_env
  ```

   Now the GAAS environment is displayed at the beginnin of each prompt line: `(GAAS)`
   To get out of the nbis environment and restore your previous environment type:

  ```bash
  deactivate
  ```

   To reactivate he GAAS environment at any time, just type:

   ```bash
   source ~/git/GAAS/profiles/activate_rackham_env
   ```

### Setup your general working place    
Then you will prepare your general working place.  

   * Move to the place where you will work  

   ```bash
   cd /proj/g2019006/nobackup/
   ```

   * Create your private place where all the magic will happen.  

   ```bash
   mkdir $USER
   cd $USER
   ```

   * get the data  
   Data needed for the exercices have to be copied locally.  

   ```bash
   cp -r /sw/courses/annotation/2019/data .
   ```
