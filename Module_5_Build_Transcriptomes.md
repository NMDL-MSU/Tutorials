
---
title: Assemble Potential Transcripts
author: Deborah Velez-Irizarry
date: Updated Jan 8 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
In this tutorial, we will cover how to assemble each bam file into potential transcripts. 

![](https://user-images.githubusercontent.com/44003875/71834718-df159400-307d-11ea-927f-8aaf83ab9482.png)

This tutorial is following the protocol presented in nature protocols:

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL (2016).
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown.](https://www.ncbi.nlm.nih.gov/pubmed/27560171)
Nature Protocols, 11(9):1650-1667 

![](https://user-images.githubusercontent.com/44003875/71835383-69aac300-307f-11ea-9da4-c3c814a36a72.png)


### Connect to HPC system
If you are using Ubuntu, open up the terminal with `Ctrl`+`Alt`+`T`. 
From the command line, log in to HPC using ssh. 

```bash
ssh -YX username@gateway.hpcc.msu.edu
```

If you are using the remote desktop environment login to your terminal 
through the Web-based remote desktop: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)

### Sort BAM Files

SAMtools has several utilities for sequence alignment and mapping for SAM/BAM format. 
It is also commonly used to call variants from RNAseq data. The [SamTools Mannual](http://www.htslib.org/doc/samtools.1.html)
can be of use if you wish to look at specific regions in your BAM files or merge files 
together. The SAM flag in `samtools view` is used to index your reads based on certain properties,
such as read first pair, read mapped in proper pair, read rever strand or any combination. 
The [Picard flags page](https://broadinstitute.github.io/picard/explain-flags.html) is a great resource to 
obtain a SAM Flag based on user specified properties. In this tutorial we will use samtools to 
calculate dapth of coverage. In this step we will use SAMtools to sort the BAM files containing the unique 
reads.


```bash
nano $HOME/RNAseq_Pipeline/SortBAM.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: SortBAM.sh
#   Directory code: $HOME/RNAseq_Pipeline/HISAT2
#   Date: January 8, 2020
#   Description: Sort unique BAM files
#   Run: bash SortBAM.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/HISAT2
#
#   Output files to directory:
#       $SCRATCH/HISAT2
#==============================================================================

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Work Directory
dir=$HOME/RNAseq_Pipeline/HISAT2

# Input Directory
inp=$SCRATCH/HISAT2

# Output Directory
out=$SCRATCH/HISAT2

# Qstat directory
mkdir $HOME/RNAseq_Pipeline/HISAT2/qstat/sort
qstat=$HOME/RNAseq_Pipeline/HISAT2/qstat/sort

# Move script to directory
mv $HOME/RNAseq_Pipeline/SortBAM.sh $dir

# Write bash script for each animal
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=50G
#SBATCH -J '${anim[$i]}_sortBAM'
#SBATCH -o '${anim[$i]}_sortBAM.o%j'

# Work Directory
cd '$out'

# Load required modules
module load SAMtools/1.9

# Module List
module list

# Sort BAM files
samtools sort -o '$out/${anim[$i]}'_uniq_sorted.bam '$inp/${anim[$i]}'_uniq.bam   

# Run statistics
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}.qsub

# Submit script to hpcc
cd $qstat
sbatch ${anim[$i]}.qsub

done
```
> Run master script.

```bash
bash $HOME/RNAseq_Pipeline/SortBAM.sh
```

Check submitted jobs `sq`. When all jobs have completed, check for errors before 
procceding to the next step.

```bash
cd $HOME/RNAseq_Pipeline/HISAT2/qstat/sort
checkJobs
```

### Obtain depth of coverage

Before builing the transcriptome, we should get an average estimate of the depth of coverage 
per animal. The average depth of coverage (focusing on covered bases) is a ratio of the total 
number of counts per covered regions over the total number of nucleotides sequenced. 
This will be estimated per bam file using Samtools. Another estimate of coverage is the 
average X coverage were instead of dividing by the total number of nucleotides sequenced, you
divide by the total genome size.

```bash
nano $HOME/RNAseq_Pipeline/depth.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: depth.sh
#   Directory code: $HOME/RNAseq_Pipeline/Depth/depth.sh
#   Date: January 8, 2020
#   Description: Obtain read depth for each animal unique mapped reads.
#   Run: bash depth.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/HISAT2
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Depth
#
#   Output files:
#       uniq_depth.txt
#		uniq_X_coverage.txt
#==============================================================================

# Input Directory
dir=$SCRATCH/HISAT2

# Output Directory
mkdir $HOME/RNAseq_Pipeline/Depth
out=$HOME/RNAseq_Pipeline/Depth

# Chromosome directory
mkdir $HOME/RNAseq_Pipeline/Depth/Chrom
chr=$HOME/RNAseq_Pipeline/Depth/Chrom

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Move bash script to Depth directory
mv $HOME/RNAseq_Pipeline/depth.sh $out

# Load SamTools Module
module load SAMtools/1.9

# Add column names to out files
echo Animal TotalCovered NucleotidesSequenced > $out/uniq_depth.txt
sed -e 's/ /\t/' $out/uniq_depth.txt > tmp; mv tmp $out/uniq_depth.txt

echo Animal TotalCovered TotalGenomeSize > $out/uniq_X_coverage.txt
sed -e 's/ /\t/' $out/uniq_X_coverage.txt > tmp; mv tmp $out/uniq_X_coverage.txt

# Calculate depth interactively per animal
for ((i=0; i<${#anim[@]} ; i++ )) do

# Work directory
cd $dir

# Generate average coverage per base sequenced per animal
samtools depth ${anim[$i]}_uniq.bam  |  awk '{sum+=$3} END { print '${anim[$i]}' "\t" sum  "\t" NR "\t" ,sum/NR}' >> $out/uniq_depth.txt

# Generate average X coverage per genome size per animal
G=(`samtools view -H ${anim[$i]}_uniq.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`)
samtools depth ${anim[$i]}_uniq.bam  |  awk '{sum+=$3} END { print '${anim[$i]}' "\t" sum  "\t" '$G' "\t" ,sum/NR}' >> $out/uniq_X_coverage.txt

# Generate depth per chromosome
samtools idxstats ${anim[$i]}_uniq.bam | awk '"'"'{print $1" "$3}'"'"' > $chr/${anim[$i]}_uniq_chr_depth.txt
echo chrAll `grep chr $chr/${anim[$i]}_uniq_chr_depth.txt | cut -f2 -d' ' | paste -s -d+ | bc` >> $chr/${anim[$i]}_uniq_chr_depth.txt
sed 's/ /\t/g' $chr/${anim[$i]}_uniq_chr_depth.txt > ${anim[$i]}; mv ${anim[$i]} $chr/${anim[$i]}_uniq_chr_depth.txt

done
```

> Run script interactively.

```bash
bash $HOME/RNAseq_Pipeline/Depth/depth.sh
```

This script will run interactively in your terminal by using a `for` loop. Once finished
review the average coverage per seuqenced base:

```bash
less $HOME/RNAseq_Pipeline/Depth/uniq_depth.txt
```

And the average X coverage

```bash
less $HOME/RNAseq_Pipeline/Depth/uniq_X_coverage.txt
```

The depth per chromosome was also calculated per animal. These text files were saved in the
`Chrom` directory.

```bash
cd $HOME/RNAseq_Pipeline/Depth/Chrom
cat `ls *_uniq_chr_depth.txt | head -1`
```

### Build poteintial transcripts

StringTie will be used to assemble transcriptomes per file. Refer to the [StringTie mannual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
for futher options that ones used in this tutorial. 

```bash
nano $HOME/RNAseq_Pipeline/StringTie.sh
```

> Copy the following script and paste in the terminal editor window.

