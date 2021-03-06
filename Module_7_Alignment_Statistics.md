---
title: Alignment Statistics RNA-seq Pipeline
author: Deborah Velez-Irizarry
date: Updated Jan 17 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
In this tutorial, we will review the performance of the RNA-seq pipeline. 
The descriptive analysis, run in R, will take as input the summary files generated 
by each step of the pipeline and evaluate the number of reads retained and dropped
across each step. 

![](https://user-images.githubusercontent.com/44003875/72625890-a5634980-3917-11ea-8831-b182dafcb1c8.png)

We will compare this pipeline with the previous pipeline using TopHat2 and Cufflinks:

![](https://user-images.githubusercontent.com/44003875/72625995-d774ab80-3917-11ea-97f6-c0a1b97cffee.png)


### Run R

The following script is a short program that will be used to run all of your
R scripts. The program create a directory using the name of your R script
and will move the script to that directory and submit it to HPCC. All the output
generated by this script will be saved in the newly created directory. 
Save this script in your home directory in a new directory called bin:

```bash
mkdir $HOME/bin
cd $HOME/bin
nano RunR
```

> Copy the following script and paste in the terminal editor window.


```bash
#!/usr/bin/env bash

#####################################################################
# Script to run knitR scripts on HPC and render html
#####################################################################

response=$@

File=`echo $response | cut -f1 -d' '`
Dir=`pwd`
Name=`echo $response | cut -f1 -d' ' | cut -f1 -d.`
Res=`echo $response | cut -f2 -d' '`
Title=`echo $response | cut -f2 -d+`
node=`echo $Res | cut -f1 -d, | cut -f2 -d=`
cpus=`echo $Res | cut -f2 -d, | cut -f2 -d=`
time=`echo $Res | cut -f3 -d, | cut -f2 -d=`
mem=`echo $Res | cut -f4 -d, | cut -f2 -d=`
D=`date`
A=`echo $USER`

mkdir $Dir/$Name
mv $Dir/$File $Dir/$Name
cd $Dir/$Name

echo '---
title: '$Title'
author: '$A'
date: '$D'
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---' > $Dir/$Name/.metadata

echo '#!/bin/bash
#SBATCH --nodes='$node'
#SBATCH --cpus-per-task='$cpus'
#SBATCH --time='$time'
#SBATCH --mem='$mem'
#SBATCH -C intel18
#SBATCH -J '$Name'
#SBATCH -o '$Name.o%j'

#=====================================================================
# This script runs: '$File'
# Submited on: '$D'
#=====================================================================

# Work Directory
cd '$Dir/$Name'

# Run R Code
'R -e "'library("'"knitr"'");knitr::spin ("'"'$File'"'", precious=TRUE)'"'

# Add metadata to md file
cat .metadata '$Name'.Rmd > tmp; mv tmp '$Name'.Rmd
cat .metadata '$Name'.md > tmp; mv tmp '$Name'.md

# Render HTML
'R -e "'library("'"prettydoc"'");rmarkdown::render ("'"'$Name.Rmd'"'")'"'

# Job Details
echo 'Job Details'
scontrol show job $SLURM_JOB_ID' > $Name.qsub

qsub $Name.qsub
```


### Descriptive Analysis

The following R script will take as input the summary files generated for each step of the RNA-seq pipeline.
For each step summary statistics will be calculated.

```bash
cd $HOME/RNAseq_Pipeline
nano Descriptive_Analysis.R
```

> Copy the following script and paste in the terminal editor window.


```bash
#' ### Description:  
#' Alignment statistics for RNA-seq Pipeline.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;$HOME/RNAseq_Pipeline  
#'  
#' Directory/File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis/Descriptive_Analysis.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Trimmomatic/trimmomatic_rst.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HISAT2/summary_alignment.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/summary_htseq.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/paired_counts.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/single_counts.txt  
#'  
#' **Output files:**  
#'  
#' Directory/File:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;$HOME/RNAseq_Pipeline/Descriptive_Analysis/retained_read_stats.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;$HOME/RNAseq_Pipeline/HTSeq/htseq_counts.txt  
#'  
#' ***  

#' ### Code
#' Clear Environment
rm(list=ls())

#' **Session Information**
sessionInfo()

#' **Input Directory**
dir <- "$HOME/RNAseq_Pipeline"
dir <- "/mnt/research/NMDL/RER_Thoroughbred/RNAseq_Pipeline"

#' ### Animal IDs
animID <- rownames(read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"),
    header=TRUE, row.names=1))[1:23]
animID

#' ### Load RNA-Seq Pipeline Statistics:
#' **Trimmomatic**
adpt <- read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"),
    header=TRUE, row.names=1)
adpt <- adpt[animID,]
dim(adpt)


#' **HISAT2**
hisat2 <- read.table(paste(dir,"HISAT2", "summary_alignment.txt", sep="/"),
    header=TRUE, row.names=1)
hisat2 <- hisat2[animID,]
dim(hisat2)


#' **HTSeq**

#' > Summary HTSeq table
sumHTSeq <- read.table(paste(dir, "HTSeq", "summary_htseq.txt", sep="/"),
    header=TRUE, row.names=1)
sumHTSeq <- sumHTSeq[animID,]
dim(sumHTSeq)

#' > Paired read counts table
paired <- read.table(paste(dir, "HTSeq", "paired_counts.txt", sep="/"),
    header=TRUE, row.names=1)
paired <- paired[,paste("X", animID, "_paired", sep="")]
dim(paired)

#' > Single read count table
single <- read.table(paste(dir, "HTSeq", "single_counts.txt", sep="/"),
    header=TRUE, row.names=1)
single <- single[,paste("X", animID, "_single", sep="")]
dim(single)

#'  ### Summary Function
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean",
    "Std.Dev.", "3rd Qu.", "Max.")]


#' ### Summary Trimmomatic: Adapter filtering
#' **Number of Input Read Paires (million reads)**
startPE <- adpt$input.read.pairs*2
names(startPE) <- rownames(adpt)
summSD((startPE/2)/1e6)

#' **Number of Input Read (million reads)**
summSD((startPE)/1e6)

#' **Number of retained paired reads after filtering adapter sequences (million)**
summSD(adpt$both.surviving*2/1e6)

#' Percent retained paired reads
mean((adpt$both.surviving*2) / startPE)*100

#' **Number of forward reads surviving without its pair (million)**
summSD(adpt$fwd.only.surviving/1e6)

#' **Number of reverse reads surviving without its pair (million)**
summSD(adpt$rev.only.surviving/1e6)

#' **Number of single reads surviving (million)**
summSD((adpt$fwd.only.surviving/1e6 + adpt$rev.only.surviving/1e6))

#' Percent single reads surviving
mean((adpt$fwd.only.surviving + adpt$rev.only.surviving) / startPE)*100

#' **Number of retained reads start after adapter trimming**
trimm.out <- rowSums(data.frame((adpt$both.surviving*2), + adpt$fwd.only.surviving +
    adpt$rev.only.surviving))
names(trimm.out) <- rownames(adpt)
summSD(trimm.out/1e6)

#' **Number of dropped reads (million)**
summSD(startPE - trimm.out)/1e6

#' Percent single reads surviving
mean((startPE - trimm.out) / startPE)*100

#' **Percent of retained reads after adapter trimming from total number of sequenced reads**
summSD(trimm.out/startPE)*100

#' **Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Trimmomatic.
#' Output should be zero.
sum(!rowSums(data.frame((adpt$both.surviving*2),
    (rowSums(data.frame(adpt$fwd.only.surviving, adpt$rev.only.surviving,
    adpt$dropped)*2)))) == startPE)


#' ### Summary HISAT2: Number of reads aligning to reference genome
#' **Data Check:** Input reads for HISAT2 should be equal to output reads for Trimmomatic
# Do the read input for HISAT2 reflect the read output for trimmomatic? Output should be zero.
start.top <- hisat2$total_input_reads
names(start.top) <- rownames(hisat2)
sum(sum(!start.top) == trimm.out[rownames(hisat2)])

#' **Number of paired reads aligning to the reference**
summSD(hisat2$paired_aligned_reads/1e6)

#' Percent paired reads
mean(hisat2$paired_aligned_reads/start.top)*100

#' **Number single reads aligning to the reference**
summSD(hisat2$unpaired_aligned_reads/1e6)

#' Paired single
mean(hisat2$unpaired_aligned_reads/start.top)*100

#' **Number of unmapped reads**
hisat2.drop <- (start.top - (hisat2$paired_aligned_reads +
    hisat2$unpaired_aligned_reads))
names(hisat2.drop) <- rownames(hisat2)
summSD(hisat2.drop/1e6)

#' Percent unmapped
mean(hisat2.drop/start.top)*100

#' **Total number of mapped reads**
hisat2.out <- (hisat2$total_aligned_reads)
names(hisat2.out) <- rownames(hisat2)
summSD(hisat2.out/1e6)

#' **Mapping percent from total input reads**
summSD(hisat2$alignment_rate)

#' Double check you get the same percent of mapped reads
mean(hisat2.out/start.top)*100

#' **Percent of reads retained (aligned) from total number of sequenced reads**
summSD(hisat2.out/startPE[names(hisat2.out)])*100



#' ### Samtools: Retain unique reads
#' **Number of uniquely alignned reads**
uniq.out <- hisat2$total_uniquely_aligned_reads
names(uniq.out) <- rownames(hisat2)
summSD(uniq.out/1e6)

#' Percent unique reads
mean(uniq.out/hisat2.out)*100

#' **Number of dropped reads**
summSD(hisat2.out - uniq.out)/1e6

#' Percent dropped
mean((hisat2.out - uniq.out)/hisat2.out)*100

#' **Mapping percent for unique reads from total input reads**
summSD(hisat2$unique_alignment_rate)

#' **Total number of reads for HTSeq from total number of sequenced reads**
summSD(uniq.out/startPE[names(uniq.out)])*100

#' **Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for hisat2.
#' Output should be zero.
sum(!rowSums(data.frame(hisat2.out, hisat2.drop)) == start.top)


#' ### Summary HTSeq:
#' **Number of reads processed**
processed.htseq <- rowSums(data.frame(sumHTSeq$Paired_Processed*2,
    sumHTSeq$Single_Processed))
names(processed.htseq) <- animID
summSD(processed.htseq/1e6)

#' **Number of reads processed with no feature**
no.feature.htseq <- rowSums(data.frame(sumHTSeq$Paired_No_Feature*2,
    sumHTSeq$Single_No_Feature))
names(no.feature.htseq) <- animID
summSD(no.feature.htseq/1e6)

#' **Percent of reads processed with no feature atribute (no counts obtained) per animal**
summSD(no.feature.htseq/processed.htseq*100)

#' **Percent of retained processed reads per animal**
summSD(100 - (no.feature.htseq/processed.htseq*100))


#' ### Summary expressed genes
#'  Merge paired and single read counts per animal
anim <- unlist(lapply(strsplit(colnames(paired), "_"),
    function(x) x[[1]][1]))
counts <- do.call(cbind, lapply(1:ncol(paired),
    function(x) rowSums(data.frame(paired[,x], single[,x]))))
colnames(counts) <- anim
rownames(counts) <- rownames(paired)
dim(counts)

#' **Remove HTSeq statistics from counts (laste five rows)**
tail(counts[,1:5])
counts <- counts[1:(nrow(counts)-5),]
tail(counts[,1:5])


#' **Total number of expressed genes (genes with at least one count)**
counts <- counts[rowSums(counts) > 0,]
nrow(counts)


#' ### Summary genes for differential expression analysis
#' **Number of gene transcripts with expression counts greater than 2 in all animals**
# Number of gene transcripts available for downstream analysis:
idx <- apply(counts,1, function(x) sum(x > 2) == ncol(counts))
counts <- counts[idx,]
nrow(counts)

#' **Save count HTSeq table**
write.table(counts, paste(dir, "HTSeq", "htseq_counts.txt", sep="/"),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


#' ### Save RNA-Seq Pipeline read summary
save(startPE, trimm.out, hisat2.out, uniq.out, processed.htseq, 
	no.feature.htseq, counts, file=paste(getwd(), 
		"retained_read_stats.Rdata", sep="/"))



#' ### Run R Script
#+ eval = FALSE
~/bin/RunR Descriptive_Analysis.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RNA-Seq Pipeline Statistics
```

> To run the RunR program:

```bash
cd $HOME/RNAseq_Pipeline
~/bin/RunR Descriptive_Analysis.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RNA-Seq Pipeline Statistics
```

Once the script finishes copy the HTML file to your desktop and open in your
browser. Use the information in the the HTML to fill out your pipeline summary:

![](https://user-images.githubusercontent.com/44003875/72632103-91254980-3923-11ea-8972-be7f3704f4fc.png)

I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.


