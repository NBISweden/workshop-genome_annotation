---
layout: default-overview
title: Get started
exercises: 15
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

  * **Connection to Uppmax**  
  Please connect yourself to Uppmax following those instruction [UPPMAX login instructions](uppmax_login).

  * **Create the folder structure**  
  Once connected you will create and move into the **annotation\_course** folder, where all the magic will happen.

  ```
  mkdir -p ~/annotation_course/check_assembly
  cd ~/annotation_course
  ```

  * **get the data**  
  Data that will be needed for the annotation course can be found THERE

  ```
  cp /proj/uppstore2019059/RESTOFPATH .
  ```

  * **install the GAAS repository through dedicated environment**  
  Here I will need you to clone and install our github repository indeed this repository contains NBIS annotation team scripts and will be used during this pratical.

  ```
  mkdir -p ~/git ; mv ~/git
  git clone https://github.com/NBISweden/GAAS.git
  cd GAAS
  make install 
  source ~/git/GAAS/profiles/activate_rackham_env
  ```
  
   Now the GAAS environment is displayed at the beginnin of each prompt line: `(GAAS)`
   To get out of the nbis environment and restore your previous environment type:

  ```
  deactivate
  ```
  
   To reactivate he GAAS environment at any time, just type:
   
   ```
   source ~/git/GAAS/profiles/activate_rackham_env
   

  * **Move into the proper folder to start the exercise**  
  Now move into the **check_assembly** folder and you are ready to start for this morning !
  
  ```
  cd ~/annotation_course/check_assembly
  ```
