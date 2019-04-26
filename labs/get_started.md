# Foreword:

We will for all exercises use data for the fruit fly, Drosophila melanogaster, as that is one of the currently best annotated organisms and there is plenty of high quality data available. However, working on eukaryotes can be time consuming. Even a small genome like Drosophila would take too long to run within the time we have for this course. Thus to be sure to perform the practicals in good conditions, we will use the smallest chromosome of the drosophila (chromosome 4) like it was a whole genome.


# Prerequisites

  * **Connection to Uppmax**  
Please connect yourself to Uppmax following those instruction[UPPMAX login instructions](uppmax_login).

  * **Create the folder structure**  
Once connected you will create and move into the **annotation\_course** folder, where all the magic will happen.
```
mkdir -p ~/annotation_course/check_assembly
cd ~/annotation_course
```

  * **get the data**  
Data that will be needed for the annotation course can be found THERE

```
ln -s /proj/uppstore2019059/RESTOFPATH . or copy?
```

  * **install the GAAS repository**  

  Here I will need you to clone our github indeed this github contain NBIS annotation team scripts and will be used during this pratical.

  ```
  git clone https://github.com/NBISweden/GAAS.git
  ```
  As we will be using the scripts libraries available in the git GAAS you need first to export the libraries :

  ```
  export PERL5LIB=$PERL5LIB:~/annotation_course/GAAS/annotation/
  ```


  * **Move into the proper folder to start the exercise**  
Now move into the **check_assembly** folder and you are ready to start for this morning !
```
cd ~/annotation_course/check_assembly
```
