---
title: Adapter and Quality Filtering
author: Deborah Velez-Irizarry
date: Updated Jan 27, 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
In this tutorial, we will cover Illumina adapter trimming and quality 
filtering. Before we start this step we will need to merge fastq files 
per animal.

![](https://user-images.githubusercontent.com/44003875/71802462-bb296280-302b-11ea-99fd-093405ec441b.png)

### Connect to HPC system
If you are using Ubuntu, open up the terminal with `Ctrl`+`Alt`+`T`. 
From the command line, log in to HPC using ssh. 

**SSH **
```bash
ssh -YX username@gateway.hpcc.msu.edu
```

If you are using the remote desktop environment login to your terminal 
through the Web-based remote desktop: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)

### Merge Fasta files
In our example we ran a second lane to increase the number of reads per 
animal so we will need to merge these extra reads for the first and second 
strand seperately. If you do not have multiple sequence runs per animal 
you can skip this step. 

First, save the master script to the working directory you created for 
the RNAseq pipeline.

```bash
nano $HOME/RNAseq_Pipeline/Merge_RNASeq_reads_per_animal.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#======================================================================
#   Script: Merge_RNASeq_reads_per_animal.sh
#   Directory Code: $HOME/RNAseq_Pipeline/Merged
#   Date: January 6, 2020
#   Description: Merge RNA-Seq reads for each animal into
#          a single fastq file.
#   Run: bash Merge_RNASeq_reads_per_animal.sh
#----------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/RAW/*fastq
#
#   Output files to scratch directory:
#       $SCRATCH/Merged/*fastq
#======================================================================

# Directory of fastq files
Raw=$SCRATCH/RAW

# Directory for merged fastq files
mkdir $HOME/RNAseq_Pipeline/Merged
mkdir $SCRATCH/Merged
Merg=$SCRATCH/Merged

# Batch script statistics directory
mkdir $HOME/RNAseq_Pipeline/Merged/qstat
qstat=$HOME/RNAseq_Pipeline/Merged/qstat

# Move script to work directory
mv $HOME/RNAseq_Pipeline/Merge_RNASeq_reads_per_animal.sh $HOME/RNAseq_Pipeline/Merged

# Animal IDs
cd $Raw
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Merge reads per animal
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH -J '${anim[$i]}'
#SBATCH -o '${anim[$i]}.o%j'

# Change to merged directory
cd '$Merg'

# Files to merge
echo "List of files to merge for '${anim[$i]}'"
ls -ltrh '$Raw'/'${anim[$i]}'_*

# Merge R1 fastq files for '${anim[$i]}'
R1=(`ls '$Raw'/'${anim[$i]}'_*R1*`)
for ((j=0; j<${#R1[@]} ; j++ )) do
cat ${R1[$j]} >> '${anim[$i]}'_R1_merged.fastq
done

# Merge R2 fastq files for '${anim[$i]}'
R2=(`ls '$Raw'/'${anim[$i]}'_*R2*`)
for ((j=0; j<${#R2[@]} ; j++ )) do
cat ${R2[$j]} >> '${anim[$i]}'_R2_merged.fastq
done

# qstat
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}.qsub

# Submit job to hpcc
cd $qstat
sbatch ${anim[$i]}.qsub

done
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Run the Merge master script with `bash`:

```bash
bash $HOME/RNAseq_Pipeline/Merge_RNASeq_reads_per_animal.sh
```

This script will create a directory in your scratch space called Merged. 
The merged fasta files will be saved to this directory.


### Trimmomatic

Trimmomatic is a flexable read trimming tool created for Illumina 
next-generation sequencing data. To learn more about Trimmomatic refer to:

> Bolger, A.M., Lohse, M., &amp; Usadel, B. (2014). [Trimmomatic: A flexible trimmer for Illumina Sequence Data.](https://www.ncbi.nlm.nih.gov/pubmed/24695404) Bioinformatics, btu170.  
> [Trimmomatic Website](http://www.usadellab.org/cms/index.php?page=trimmomatic)  


The parameters used to performe this step are as follows:

**ILLUMINACLIP:**  
Cuts the Illumina-sequence specific adapters. We will trim the TruSeq3-PE 
adapter sequences shown below.

> PrefixPE/1 TACACTCTTTCCCTACACGACGCTCTTCCGATCT  
> PrefixPE/2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  
> PE1 TACACTCTTTCCCTACACGACGCTCTTCCGATCT  
> PE1_rc AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA  
> PE2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  
> PE2_rc AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  

**HEADCROP**  
Cut a specific number of bases from the start of the read. We will trim 
the first 6 bases.

**LEADING**  
Cut the leading bases only if below a quality of 3

**TRAILING**  
Cut the trailing bases only if below a quality of 3

**SLIDINGWINDOW**  
Sliding window to scan read, we will use a 4-base window and trim if the 
average quality is below 15. 

**MINLEN**  
Drop reads with a minumin length of 75 bases.

  
First, save the master script to the working directory you created for 
the RNAseq pipeline tutorial.

```bash
nano $HOME/RNAseq_Pipeline/adapter_trimming.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#======================================================================
#   File: adapter_trimming.sh
#   Directory code: $HOME/RNAseq_Pipeline/Trimmomatic
#   Date: January 6, 2020
#   Description: Trim adapter sequences using Trimmomatic software
#   Run: bash adapter_trimming.sh
#----------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/Merged/*fastq
#
#   Output files to directories in scratch:
#       $SCRATCH/Trimmomatic/PE
#       $SCRATCH/Trimmomatic/SE
#======================================================================

### Raw RNA-seq Directory
raw=$SCRATCH/Merged

### Animal IDs
cd $raw
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

### Output directories
# Trimmomatic directory
mkdir $HOME/RNAseq_Pipeline/Trimmomatic
mkdir $SCRATCH/Trimmomatic
trim=$SCRATCH/Trimmomatic

# Bash qstatistics directory
mkdir $HOME/RNAseq_Pipeline/Trimmomatic/qstat
qstat=$HOME/RNAseq_Pipeline/Trimmomatic/qstat

# Paired end reads directory
mkdir $SCRATCH/Trimmomatic/PE
pe=$SCRATCH/Trimmomatic/PE

# Single end reads directory
mkdir $SCRATCH/Trimmomatic/SE
se=$SCRATCH/Trimmomatic/SE

## Move bash script to work directory
mv $HOME/RNAseq_Pipeline/adapter_trimming.sh $HOME/RNAseq_Pipeline/Trimmomatic


### Generate bash scripts to run trimmomatic per animal
for((i=0;$i<${#anim[@]};i++))
do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH -J '${anim[$i]}'
#SBATCH -o '${anim[$i]}_%j'

module load Trimmomatic

# Work directory
cd '$trim'

# Filter adaptor sequences from fastq
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE '$raw'/'${anim[$i]}'_R1_merged.fastq '$raw'/'${anim[$i]}'_R2_merged.fastq '${anim[$i]}'_R1_merged_PE.fastq '${anim[$i]}'_R1_merged_SE.fastq '${anim[$i]}'_R2_merged_PE.fastq '${anim[$i]}'_R2_merged_SE.fastq ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:6 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

# Move paired end files to PE and singles to SE
mv '${anim[$i]}'_R1_merged_PE.fastq '${anim[$i]}'_R2_merged_PE.fastq '$pe'
mv '${anim[$i]}'_R1_merged_SE.fastq '${anim[$i]}'_R2_merged_SE.fastq '$se'

# qstat
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}.qsub

# Submit job to hpcc
cd $qstat
sbatch ${anim[$i]}.qsub

done
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Submit the master script for Trimmomatic to the SLURM scheduler.

```bash
bash $HOME/RNAseq_Pipeline/adapter_trimming.sh
```

This step should take about an hour to run. Wait until all jobs have 
finished before continuing to the step. To check the jobs status use `qs`. 
When all jobs have finished run the checkjobs within the `qstat` directory 
to make sure the jobs ran with no errors. 

```bash
cd $HOME/RNAseq_Pipeline/Trimmomatic/qstat
checkJobs
```

### Generate summary for quality trimming

To keep a complete documentation of retained and dropped reads we will 
create a summary file for each step in the RNA-seq bioinformatic pipeline. 
We will extract this information from the job outfiles (one per animal or 
animal-timepoint depending on the experimental design). Let us look at how 
this information is provided by Trimmomatic:

```bash
cd $HOME/RNAseq_Pipeline/Trimmomatic/qstat
cat `ls | grep -v qsub | head -1`
```

Observe the number provided for `Input Read Pair`, `Both Surviving`, 
`Forward Only Surviving`, `Reverse Only Surviving`, and `Dropped`. 
Notice that this quality step generates both paired and single 
reads. This happens when the trimming results in a complete loss of one 
of the read paires. We will extract this information for each of the jobs 
submitted and save it to a tab delimited matrix that we will use at the 
end of the pipeline to asses retained and loss reads. 


```bash
nano $HOME/RNAseq_Pipeline/merge_rst_trim_adapters.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#===================================================================
#   File: merge_rst_trim_adapters.sh
#   Directory: $HOME/RNAseq_Pipeline/Trimmomatic
#   Date: January 6, 2020
#   Description: Check that all files were properly trimmed and merge
#                all the results.
#   Run: bash merge_rst_trim_adapters.sh
#-------------------------------------------------------------------
#   Input files in directory:
#       $HOME/RNAseq_Pipeline/Trimmomatic/qstat
#
#   Output file to directory:
#       $HOME/RNAseq_Pipeline/Trimmomatic
#
#   Output file:
#       rst_trim_adapt.txt
#       trimmomatic_rst.txt
#====================================================================

### Raw RNA-seq Directory
raw=$SCRATCH/RAW

### Animal IDs
cd $raw
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

### Output directory
trim=/$SCRATCH/Trimmomatic

### Input directory
rst=$HOME/RNAseq_Pipeline/Trimmomatic

### Qstat directory
qstat=$HOME/RNAseq_Pipeline/Trimmomatic/qstat

## Move bash script to Trimmomatic directory
mv $HOME/RNAseq_Pipeline/merge_rst_trim_adapters.sh $HOME/RNAseq_Pipeline/Trimmomatic

# Check if you have an output for each file
pe=$trim/PE
se=$trim/SE

# Paired
# If all jobs generated output you should not get any output.
for ((i=0; i<${#anim[@]} ; i++ )) do
ls $pe/${anim[$i]}_R1_merged_PE.fastq > $pe/temp
ls $pe/${anim[$i]}_R2_merged_PE.fastq > $pe/temp
done; rm $pe/temp

# Single
# If all jobs generated output you should not get any output.
for ((i=0; i<${#anim[@]} ; i++ )) do
ls $se/${anim[$i]}_R1_merged_SE* > $se/temp
ls $se/${anim[$i]}_R2_merged_SE* > $se/temp
done; rm $se/temp

# Check if all files completed successfully
# If all jobs generated output you should not get any output.
cd $qstat
for ((i=0; i<${#anim[@]} ; i++ )) do
echo `echo -n ''${anim[$i]}' ' ;  grep 'TrimmomaticPE: C' ${anim[$i]}_*` >> ../completed.txt
done

# Merge Results and save to text file
grep 'Input Read Pairs' * > $rst/rst_trim_adapt.txt

# Creat text file with trimmomatic counts for each animal
D=$rst/rst_trim_adapt.txt
cat $D | cut -f1 -d_ > id
cat $D | cut -f4 -d' ' > input.read.pairs
cat $D | cut -f7 -d' ' > both.surviving
cat $D | cut -f12 -d' ' > fwd.only.surviving
cat $D | cut -f17 -d' ' > rev.only.surviving
cat $D | cut -f20 -d' ' > dropped

echo id input.read.pairs both.surviving fwd.only.surviving rev.only.surviving dropped > trimmomatic_rst.txt
cat trimmomatic_rst.txt | tr ' ' '\t' > tmp; mv tmp $rst/trimmomatic_rst.txt
paste id input.read.pairs both.surviving fwd.only.surviving rev.only.surviving dropped >> $rst/trimmomatic_rst.txt
rm id input.read.pairs both.surviving fwd.only.surviving rev.only.surviving dropped
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Run the script using bash. 

```bash
bash $HOME/RNAseq_Pipeline/merge_rst_trim_adapters.sh
```

This bash script will be run interactively; meaning that it will not be submitted to the 
scheduler. An error generated from this script will point you to a job that did not generate 
output. If you get an error rerun the qsub script for that file. Refer to the qstat job output 
file for further instructions. If the job did not finish in the time allocated or used up more 
memory than requested you will need to edit the qsub file and increase the requested resorces 
before resubmitting the job to SLURM. 

Let us review the summary files generated.

```bash
cd $HOME/RNAseq_Pipeline/Trimmomatic
cat trimmomatic_rst.txt
```

We would expect a high percent of retained paired reads and about less than 2.5% dropped reads 
in this step. Let us review the percent of retained paired reads and dropped reads.

```bash
cd $HOME/RNAseq_Pipeline/Trimmomatic
cat rst_trim_adapt.txt | cut -f1,8,21 -d' '
```


### Assessing Filtered Read Quality: FASTQC 

We have successfully removed Illumina adapter sequences from our RNA-seq data and filtered out 
low quality bases. Let us review the quality of these filtered reads. We need to first run
`FASTQC` followed by `MultiQC` to concatenate the FASTQC summaries. 

```bash
nano $HOME/RNAseq_Pipeline/FastQC_quality_filtered.sh
```

> Copy the following script and paste in the terminal editor window.


```bash
#==============================================================================
#   File: FastQC_quality_filtered.sh
#   Directory code: $HOME/RNAseq_Pipeline/Quality/Filtered
#   Date: January 6, 2020
#   Description: Generate FastQC reports for raw sequence reads.
#   Run: bash FastQC_Filtered.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/Trimmomatic
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Quality/Filtered
#==============================================================================

### Directory input files
filter=$SCRATCH/Trimmomatic

### Directory output files
## Create output directory
mkdir $HOME/RNAseq_Pipeline/Quality/Filtered

# Output directory
out=$HOME/RNAseq_Pipeline/Quality/Filtered

## Create qstat directory for SLURM input/output files
mkdir $HOME/RNAseq_Pipeline/Quality/Filtered/qstat

# SLURM input/output directory
qstat=$HOME/RNAseq_Pipeline/Quality/Filtered/qstat

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Write bash script for each animal
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH -J '${anim[$i]}_FastQC'
#SBATCH -o '${anim[$i]}_FastQC.o%j'

module load FastQC

# Change to filtered fastq directory
cd '$filter'

# Generate quality reports for paired reads
fastqc --outdir='$out' ./PE/'${anim[$i]}'*

# Generate quality reports for paired reads
fastqc --outdir='$out' ./SE/'${anim[$i]}'*

# qstat
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}_FastQC.qsub

# Submit job to hpcc
cd $qstat
sbatch ${anim[$i]}_FastQC.qsub

done
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Run master script:

```bash
bash $HOME/RNAseq_Pipeline/FastQC_quality_filtered.sh
```

The FASTQC master will create a new directory in your Quality folder called `Filtered`. 
Let us take a look at one of the submitted scripts:

```bash
cd $HOME/RNAseq_Pipeline/Quality/Filtered/qstat
cat `ls *.qsub* | tail -1`
```

Check the status of the submitted jobs by using the show jobs shortcut:

```bash
sq
```

For now, we wait for the jobs to finish. Only when the `sq` displays zero 
running jobs can we proceed to the next step.


### Generate MultiQC Report: Filtered

Before checking the quality summaries for the raw sequence, files let us 
first make sure all our jobs ran without errors or warnings. 
Change to the `qstat` directory and run `checkJobs`.

```bash
cd $HOME/RNAseq_Pipeline/Quality/Filtered/qstat
checkJobs
```

If the jobs finished with no errors we can proceed to generate the 
MultiQC report:

```bash
nano $HOME/RNAseq_Pipeline/multiQC_quality_filtered.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH -J 'MultiQC'
#SBATCH -o '$HOME/RNAseq_Pipeline/Quality/Filtered/MultiQC/MultiQC.o%j'

#==============================================================================
#   File: multiQC_quality_filtered.sh
#   Directory code: $HOME/RNAseq_Pipeline/Quality/Filtered/MultiQC
#   Date: January 6, 2020
#   Description: Concatenate FASTQC quality reports into a single file.
#   Run: sbatch multiQC_quality_filtered.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $HOME/RNAseq_Pipeline/Quality/Filtered
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Quality/Filtered/MultiQC
#==============================================================================

### Directory input files
mkdir $HOME/RNAseq_Pipeline/Quality/Filtered/output
Qrep=$HOME/RNAseq_Pipeline/Quality/Filtered/output
mv $HOME/RNAseq_Pipeline/Quality/Filtered/*.zip $Qrep

### Directory output files
## Create MultiQC output directory
mkdir $HOME/RNAseq_Pipeline/Quality/Filtered/MultiQC

# Output directory
QC=$HOME/RNAseq_Pipeline/Quality/Filtered/MultiQC

#' Activate virtual environment
cd $HOME/RNAseq_Pipeline/Quality/Raw/MultiQC
source QC/bin/activate

# Generate a summary file
cd $QC
multiqc $Qrep
deactivate

# Move bash script to work directory
mv $HOME/RNAseq_Pipeline/multiQC_quality_filtered.sh $QC

# Job details
scontrol show job $SLURM_JOB_ID
```

Run MultiQC on the filtered FASTQ files.

```bash
sbatch $HOME/RNAseq_pipeline/multiQC_raw.sh
```

Download the MultiQC report for quality filtered reads and compare with the previous MultiQC
report generated for the raw sequence reads. You should see an improvement in the base quality
phred scores and see no adpater sequences.

I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.


