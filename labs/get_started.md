# Foreword:

We will for all exercises use data for the fruit fly, Drosophila melanogaster, as that is one of the currently best annotated organisms and there is plenty of high quality data available. However, working on eukaryotes can be time consuming. Even a small genome like Drosophila would take too long to run within the time we have for this course. Thus to be sure to perform the practicals in good conditions, we will use the smallest chromosome of the drosophila (chromosome 4) like it was a whole genome.


# Prerequisites

FOLLOWING SHOULD BE CHANGED with connection to uppmax, cp data, get GAAS? add to PATH? module load busco and augustus, what about swissprot/ uniprot? other


  * **Connection to your virtual machine**  
Before going into the exercises below you need to connect to your virtual machine Ubuntu 16.04 following the instruction we will provide you.

  * **Create the folder structure**  
Once connected you will create and move into the **annotation\_course** folder, where all the magic will happen.
```
mkdir -p ~/annotation_course/practical1
cd ~/annotation_course
```

  * **List of tools needed. For your convenience they all have been pre-installed.**  

    * BUSCO
    * augustsus
    * GAAS repository

  * **Download the data**  
You must download the archive of the data and uncompress it (it could take few minutes).
```
wget https://u-ip-81-109.hpc2n.umu.se/tickets/7mIStX-Y-zjj_XPzI-iYQni2_0LVBSdBtHf_vhiA_Zk/data.tar.gz/download
tar xzvf download
rm download
```

  * **Move into the proper folder to start the excercice**  
Now move into the **practical1** folder and you are ready to start for this morning !
```
cd ~/annotation_course/practical1
```

## Install GAAS