---
title: Quantify Gene Expression
author: Deborah Velez-Irizarry
date: Updated Jan 14 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
In this tutorial, we will quantify individual gene profiles using the 
sorted BAM files and merged GTF. 

![](https://user-images.githubusercontent.com/44003875/72337184-f9033800-3676-11ea-8ec7-932db39a3631.png)


### HTSeq

For each file input we will obtain a count of all the reads aligning to a gene
using the annotation features (merged GTF) for each gene. The hene counts are 
generated for paired and single reads seperately. There are three modes in HTSeq,
refer to the mannual linked below, we will use `intersection-nonempty`:


![](https://user-images.githubusercontent.com/44003875/72337262-26e87c80-3677-11ea-951e-d2d0714cd715.png)

The options used in this tutorial for HTSeq are:
> `-m` specify mode intersection-nonempty  
> `-s reverse` specify the data is strand specific with the first 
reads on the opposit strand.  

Reference paper:

Simon Anders, Paul Theodor Pyl, Wolfgang Huber
HTSeq — [A Python framework to work with high-throughput sequencing data](https://academic.oup.com/bioinformatics/article/31/2/166/2366196)
Bioinformatics (2014)

Mannual:

Refer to the [mannual](https://htseq.readthedocs.io/en/release_0.11.1/count.html) 
for futher details on HTSeq modes to count features.

```bash
nano $HOME/RNAseq_Pipeline/htseq_counts_merged.sh
```

> Copy the following script and paste in the terminal editor window.


```bash
#========================================================================
#   File: htseq_counts_merged.sh
#   Directory code: $HOME/RNAseq_Pipeline/htseq_counts_merged.sh
#   Date: Jan 14 2020
#   Description: Count the number of reads in each samples (uniq.bam)
#           that map back to the merged.gtf
#       Run: bash htseq_counts_merged.sh
#-----------------------------------------------------------------------
#       Input files:
#           $SCRATCH/HISAT2/*_uniq.bam
#           $HOME/RNAseq_Pipeline/StringTie/Merged/merged.gtf
#
#       Output file to directory:
#           $HOME/RNAseq_Pipeline/HTSeq/Counts
#
#       Output file:
#           *_counts.txt
#=======================================================================

# Input Directories
bam=$SCRATCH/HISAT2
gtf=$HOME/RNAseq_Pipeline/StringTie/Merged/merged.gtf

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Work Directory
mkdir $HOME/RNAseq_Pipeline/HTSeq
home=$HOME/RNAseq_Pipeline/HTSeq

# Output directory
mkdir $HOME/RNAseq_Pipeline/HTSeq/Counts
mkdir $HOME/RNAseq_Pipeline/HTSeq/Counts/qstat
out=$HOME/RNAseq_Pipeline/HTSeq/Counts
qstst=$HOME/RNAseq_Pipeline/HTSeq/Counts/qstat

# Create virtual environment to install the required python packages
virtualenv $home/htseq
source $home/htseq/bin/activate
pip install numpy
pip install HTSeq
deactivate
act=$home/htseq/bin/activate

# Move bash script to work directory
mv $HOME/RNAseq_Pipeline/htseq_counts_merged.sh $home

# Generate Gene Counts with HTSeq
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=50G
#SBATCH -J '${anim[$i]}_htseq'
#SBATCH -o '${anim[$i]}_htseq.o%j'

# Load Modules
module load SAMtools/1.9
module list

# Divide bam file for '${anim[$i]}' into paired and single reads
TMP=/mnt/local/$SLURM_JOB_ID
samtools view -bf 1 '$bam'/'${anim[$i]}'_uniq.bam > $TMP/'${anim[$i]}'_uniq_paired.bam
samtools view -bF 1 '$bam'/'${anim[$i]}'_uniq.bam > $TMP/'${anim[$i]}'_uniq_single.bam

# Sort BAM files and convirt to SAM
typ=(`ls $TMP/*.bam | cut -f3 -d_ | cut -f1 -d.`)
for ((i=0; i<${#typ[@]} ; i++ )) do
samtools sort -n $TMP/'${anim[$i]}'_uniq_${typ[$i]}.bam -o $TMP/'${anim[$i]}'_uniq_${typ[$i]}_sorted.bam
samtools view -h $TMP/'${anim[$i]}'_uniq_${typ[$i]}_sorted.bam > $TMP/'${anim[$i]}'_uniq_${typ[$i]}_sorted.sam
done

# Run HTSeq
source '$act'

# Paired reads
htseq-count -m intersection-nonempty -s reverse $TMP/'${anim[$i]}'_uniq_paired_sorted.sam \
    '$gtf' \
    > '$out'/'${anim[$i]}'_paired_counts.txt

# Single reads
htseq-count -m intersection-nonempty -s reverse $TMP/'${anim[$i]}'_uniq_single_sorted.sam \
    '$gtf' \
    > '$out'/'${anim[$i]}'_single_counts.txt

deactivate

# Job details
echo 'Job Details'
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}_htseq.qsub

# Submit script to hpcc
cd $qstat
sbatch ${anim[$i]}_htseq.qsub

done
```

> Run master script.

```bash
bash $HOME/RNAseq_Pipeline/htseq_counts_merged.sh
```

Check submitted jobs `sq`. When all jobs have completed, 
check for errors before procceding to the next step.

```bash
cd $HOME/RNAseq_Pipeline/HTSeq/Counts/qstat
checkJobs
```


### Concatenate Counts Files

The HTSeq script creates a single file per animal. The following script
will take ach file as input and output a single file with counts
for all animal.

```bash
nano $HOME/RNAseq_Pipeline/merge_counts.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#========================================================================
#   File: $HOME/RNAseq_Pipeline/merge_counts.sh
#   Date: Jan 14 2020
#   Description: Save HTSeq summary: number of processed reads and number of
#          reads with no feature. Merge counts per animal for paired reads
#          into one file and single read counts into another.
#   Run: bash merge_counts.sh
#------------------------------------------------------------------------
#       Input files in directory:
#           $HOME/RNAseq_Pipeline/HTSeq
#
#       Input Folders:
#           Counts
#           qstat
#
#       Output file to directory:
#           $HOME/RNAseq_Pipeline/HTSeq
#
#      Output File:
#           summary_htseq.txt
#           paired_counts.txt
#           single_counts.txt
#========================================================================


# Input Directory
inp=$HOME/RNAseq_Pipeline/HTSeq/Counts

# Output Directory
out=$HOME/RNAseq_Pipeline/HTSeq

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Save number of procesed reads and number of reads with no features per animal
cd $inp/qstat
echo Animal Paired_Processed Paired_No_Feature Single_Processed Single_No_Feature > $out/summary_htseq.txt
echo ""  >> $out/summary_htseq.txt
for ((i=0; i<${#anim[@]} ; i++ )) do
    NF=`grep feature $inp/${anim[$i]}_*_counts.txt | cut -f2`
    P=`grep processed ${anim[$i]}_htseq.o* | grep -v 00000 | grep SAM | cut -f1 -d' '`
    echo `echo ${anim[$i]}` `echo $P | cut -f1 -d' '` `echo $NF | cut -f1 -d' '` \
    `echo $P | cut -f2 -d' '` `echo $NF | cut -f2 -d' '` >> $out/summary_htseq.txt
done
sed 's/ /\t/g' $out/summary_htseq.txt > tmp; mv tmp $out/summary_htseq.txt

# Merge paired counts per animal into a single file
cd $inp
cat ${anim[1]}_paired_counts.txt | cut -f1 > paired
echo XlocID > pairedID
for ((i=0; i<${#anim[[@]} ; i++ )) do
    ls ${anim[$i]}_paired_counts.txt | cut -f1,2 -d_  >> pairedID
    cat ${anim[$i]}_paired_counts.txt | cut -f2 > add
    paste paired add > tmp; mv tmp paired
    rm add
done
cat pairedID | tr "\n" " " > $out/paired_counts.txt
echo "" >> $out/paired_counts.txt
cat paired >> $out/paired_counts.txt
sed 's/ /\t/g' $out/paired_counts.txt > tmp; mv tmp $out/paired_counts.txt

# Merge single counts per animal into a single file
cat ${anim[1]}_single_counts.txt | cut -f1 > single
echo XlocID > singleID
for ((i=0; i<${#anim[@]} ; i++ )) do
    ls ${anim[$i]}_single_counts.txt | cut -f1,2 -d_ >> singleID
    cat ${anim[$i]}_single_counts.txt | cut -f2 > add
    paste single add > tmp; mv tmp single
    rm add
done
cat singleID | tr "\n" " " > $out/single_counts.txt
echo "" >> $out/single_counts.txt
cat single >> $out/single_counts.txt
sed 's/ /\t/g' $out/single_counts.txt > tmp; mv tmp $out/single_counts.txt

# Clear tmp files
rm animID.txt pairedID singleID paired single
```

> Run script interactively and review output files.

```bash
bash $HOME/RNAseq_Pipeline/merge_counts.sh
```

I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.


