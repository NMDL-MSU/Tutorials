---
title: Assessing Read Quality
author: Deborah Velez-Irizarry
date: Updated Jan 6 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
The first tutorial reviewed basic LINUX/UNIX commands and how to navigate HPCC. 
In this tutorial, we will go over the first step in an RNA-Seq pipeline, assessing read quality.

![](https://user-images.githubusercontent.com/44003875/71802323-6ab20500-302b-11ea-93e5-6c8eabab5caf.png)

### Connect to HPCC system
If you are using a Mac or Linux machine, open up the terminal and 
from the command line, log in to HPCC using ssh. 

```bash
ssh -YX username@gateway.hpcc.msu.edu
```

If you are using the remote desktop environment login to your terminal 
through the Web-based remote desktop: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)


### Create the RNA-Seq pipeline directory
The first thing you want to do when starting a new project is creating a working directory. 
Let us create a folder in your home space for the RNA-Seq pipeline:

```bash
mkdir $HOME/RNAseq_Pipeline
```

Our research space contains the sequence files for this tutorial.

```bash
cd /mnt/research/NMDL/Tutorials/RAW/MFM
```

There are a total of 48 fastq files in the RAW directory for 12 horses. 
Check the total number of files by first listing `ls` all the files in the directory and 
passing that list `|` to the word count `wc` command to count the number of lines in the 
list (lines `-l`):

```bash
ls | wc -l
```

We can obtain the total number of horses with sequence files by extracting the horse IDs 
from the file names `cut -f1 -d_`, retaining only the unique `uniq` IDs, and counting 
these unique names:

```bash
ls | cut -f1 -d_ | uniq | wc -l
```

To get the animal IDs for the sequence files, we can use the same command as before except 
we do not pipe `|` the unique IDs to the word count command. We will be using this command 
frequently when running `for` loops.

```bash
ls | cut -f1 -d_ | uniq
```

Our recent RNA-Seq experiment yielded between 30-50 million paired-end reads for a total 
of 316 fastq files and 450GB of storage space (compressed). To avoid the overutilization 
of your home directory space, we are going to copy the fastq files for this tutorial to 
your scratch directory. The sequence files copied to scratch are compressed. Before moving 
on to the next step, we need to decompress each file. Keep in mind each sequence file can be 
between 1 to 19 GB (decompressed). Coping and decompressing each file can take between 10 
minutes to one hour. We do not want to wait several hours for this simple step to finish. 
It takes a total of 9 hours to decompress the 48 files and 80 hours for all 316 files. 
Of course, I did not wait that long. The beauty of a high-performance computer is the 
ability to run multiple jobs simultaneously. It took less than 1.5 hours to run all the jobs. 
There are several ways to parallelize jobs; the simplest is to write a master script that will 
generate several smaller scripts, one for each file, and submit it to SLURM, HPCC scheduling system. 
First, save the master script to the working directory you created for the RNAseq pipeline tutorial.

```bash
nano $HOME/RNAseq_Pipeline/copy_fastq_scratch.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: copy_fastq_scratch.sh
#   Directory code: $HOME/RNAseq_Pipeline
#   Date: Janurary 6, 2020
#   Description: Copy fastq files to scratch and unzip each file.
#   Run: bash copy_fastq_scratch.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       /mnt/research/NMDL/Tutorials/RAW/*.fastq
#
#   Output files to directory:
#       $SCRATCH/RAW/*.fastq
#       $SCRATCH/RAW/qstat
#==============================================================================

### Directory input files
raw=/mnt/research/NMDL/Tutorials/RAW/MFM


### Directory output files
## Create RAW folder in scratch
mkdir $SCRATCH/RAW

# Output directory
out=$SCRATCH/RAW

## Create qstat directory for SLURM input/output files
mkdir $SCRATCH/RAW/qstat

# SLURM input/output directory
qstat=$SCRATCH/RAW/qstat


### Fastq Files
cd $raw
fl=(`ls *.fastq.gz`)

### File IDs
Id=(`ls *.fastq.gz | cut -f1,2,3,4 -d_`)

# Write bash script for each file and submit to SLURM
for ((i=0; i<${#fl[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=5G
#SBATCH -J '${Id[$i]}'
#SBATCH -o '${Id[$i]}.o%j'

# Copy fastq files to scratch
scp -p '$raw'/'${fl[$i]}' '$out'

# Gunzip fastq file
gunzip '$out'/'${fl[$i]}'

# Write job details to output file
scontrol show job $SLURM_JOB_ID' > $qstat/${Id[$i]}.qsub

# Submit job to hpcc
cd $qstat
sbatch ${Id[$i]}.qsub

done
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

This bash script will first create a few directories in your scratch space. One will contain 
the raw sequence files `RAW` and the another will save the job scripts and job details `qstat`. 
These job scripts are important to save in case an error occurs with one of your jobs, 
in which case you can review the error in the job detail, troubleshoot, and resubmit the job. 
We will review ways to check multiple job scripts for errors later in this tutorial. After 
creating the directories, a `for` loop is run on the fastq file names to create a single script 
per file. The `qsub` created in this step is saved in the `qstat` directory and is then 
submitted `sbatch` to the scheduler. The qsub script contains three commands:

> Secure copy `scp` fastq files to scratch directory.
> Decompress the fastq file.
> Save job details.

The secure copy command `scp` has the `-p` option to preserve modification times from the 
original file. To use `scp` specify the file of interest before the directory that will house 
the new copy. For example:

```bash
scp -p /mnt/research/NMDL/file.txt $HOME
```

The fastq files are compressed using gzip and contain the suffix `.gz`. 
To decompress the sequence file, we use `gunzip`. For example:

```bash
gunzip sequence_file.fastq.gz
```

Using the `gzip` command, you can compress a file. For example:

```bash
gzip sequence_file.fastq
```

We reviewed what the master script does. Now we can run it on the command line using `bash`.

```bash
bash copy_fastq_scratch.sh
```

You will see `Submitted batch job X` printed to the command line, one for each script 
summited to SLURM (48 total). Once all scripts are submitted you can check the status of 
submitted jobs by using the following command, replacing username with netid:

```bash
squeue -u username
```

**NOTE:** One can avoid having to type `squeue -u username` to check jobs by using `powertools` 
that have shortcuts to many highly used commands. To use the power tools module copy the 
following command to the terminal:

```bash
module load powertools
```

Now you can check the status of jobs by using the following shortcut:

```bash
sq
```

Adding powertools to your bashrc file in your home directory will ensure this command
is run everytime you log in. There are additional shortcuts we can add to this script
for ease of use. 

```bash
nano $HOME/.bashrc
```

> Copy the following lines and paste in the terminal editor window.

```bash
# Environment
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load R/3.5.1-X11-20180604
module load powertools
module load Pandoc
module load UDUNITS/2.2.26

# Permissions for group
umask g+rw

# Shortcut to ls
alias l='ls -ltrh'

# Directory space
alias space='du -shc ./*/'

# ssh to dev-intel18
alias log='ssh dev-intel18'

# Check submitted jobs and save resource use to qstat
alias checkJobs='bash $HOME/.checkJobs'
```

> Save the .bashrc file by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Throughout these tutorials we will be submitting several jobs to the SLURM schedular. 
I have created a script `.checkJobs` to review multiple output files from SLURM for errors or warnings. 
To use this script you will need to save it to your home directory.

```bash
nano $HOME/.checkJobs
```

> Copy the following lines and paste in the terminal editor window.

```bash
#======================================================================================
# Script to check jobs submitted to SLURM and save job resource use to output file.
# NOTE: output/error file for jobs run with SLURM should be saved in qstat folder
#======================================================================================

# Job qstat output directory: will print name of parent folder
dir=(`pwd`)
echo ''
echo Checking jobs for: `echo $dir | sed 's=/= =g' | tr ' ' '\n' | tail -2 | head -1`
echo ''

# qstat files
qstat=(`ls *.o*`)

# Job Ids
job=(`ls *.o* | cut -f2 -d. | sed 's/o//'`)

# Get state of job
for ((i=0; i<${#qstat[@]} ; i++ )) do
    echo ${qstat[$i]} `sacct -j ${job[$i]} --format=State | sed '$!d' | sed 's/   //'` >> .jobState
done

# Check for errors, warnings and scripts that are still running
echo ''
if grep RUNNING .jobState ; then
    echo ''
    echo These jobs are still running.
    echo ''

elif grep -i 'ERROR' * ; then
  echo ''
  echo There were problems in these scripts.
  echo ''

elif grep -i 'WARN' * ; then
  echo ''
  echo There were minor problems in these scripts.
  echo ''

else
  echo All jobs completed without errors.
  echo ''
fi

# Check if job output already contains run summary
done=(`grep SlurmJob *.o* | cut -f1 -d:`)

# Remove files with run status
for ((i=0; i<${#done[@]} ; i++ )) do
    sed -i '/'${done[$i]}'/d' .jobState
done

# Completed qstat files
qstat=(`grep COMPLETED .jobState | cut -f1 -d' '`)

# Job Ids
job=(`grep COMPLETED .jobState | cut -f1 -d' ' | cut -f2 -d. | sed 's/o//'`)

# Add run summary to job output file
for ((i=0; i<${#qstat[@]} ; i++ )) do
    echo '' >> ${qstat[$i]}
    echo SlurmJob Resource Use >> ${qstat[$i]}
    sacct -j ${job[$i]} --format=JobID,Jobname,Start,End,Elapsed,MaxRSS,State  >> ${qstat[$i]}
done

# Remove temporary files
rm .jobState
```
> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.


**NOTE:** Wait till all jobs have completed before continuing to the next step.


### FASTQ files and Phred Quality Scores

RNA-seq data are saved as short reads using the fastq file formate that has four lines per 
sequenced read. Let us review the first eight lines of a fastq file.

```bash
cd $SCRATCH/RAW
head -8 `ls *.fastq | head -1`
```

The first line starts with @ followed by a sequence identifier:

| Text | Description |
| --- | --- |
| `K00392` | instrument name |
| `129` | run ID |
| `HW7HLBBXX` | flowcell ID |
| `4` | flowcell lane |
| `1101` | tile number within the flowcell lane |
| `13859` | x coordinate of the cluster within the tile |
| `1349` | y coordinate of the cluster within the tile |
| `1` | member of a pair |
| `N` | for filtered reads Y, N otherwise |
| `0` | none of the control bits are on; otherwise, it is an even number |
| `GAATTCGC+TCAGAGCC` | index sequence |


The second line is the read sequence, the third has the `+` symbol, and the fourth contains 
the quality values for the read sequence. Each value inline four corresponds to a 
SANGER Phred +33 quality score.

![](https://user-images.githubusercontent.com/44003875/48494470-f2233000-e7fb-11e8-823c-02fc8c9cb939.png)


**For example:** The following table shows 12 bases of a read, the quality value for each position, and the resulting PHRED quality score.

| Line Position | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  |
| ---      | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Sequence** | A   | C   | A   | C   | G   | T   | C   | T   | G   | A   | A   | C   |
| **Quality**  | J   | <   | 7   | F   | A   | F   | A   | )   | 7   | 7   | -   | ?   |
| **PHRED**    | 41  | 27  | 22  | 37  | 32  | 37  | 32  | 8   | 22  | 22  | 12  | 30  |


By now, you are probably questioning what a Phred quality score is and how it is generated per 
base sequenced. A Phred quality score is the base calling error probabilities. For example, 
a quality score of 40 equals a 1 in 10,000 probability or 99.99% accuracy. To calculate the 
Phred quality score, the core facility at MSU incorporates PhiX control libraries during 
sequencing. The sequence generated is compared to the control library to estimate the error 
probabilities. PhiX has a defined viral genome sequence with a balanced representation 
of A, T, G, and C nucleotide and small genome size. These characteristics enable PhiX libraries 
to be used for alignment and estimation of error rates as quality control for cluster generation, 
sequencing, and alignment. The Real-Time Analysis [(RTA)](https://support.illumina.com/sequencing/sequencing_software/real-time_analysis_rta.html) 
software is used to perform base calls and quality scoring. For more information, 
see Illumina [Quality Scores](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) and [PhiX](https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) 
technical notes.


### Assessing Read Quality: FASTQC

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) software provides summary 
graphs to asses the quality of your data, including per base assessment of quality scores, 
base sequence, and GC content, and provides a review of sequence length distribution, 
duplication levels, and adapter content. We do not need to worry about installing the FASTQC, 
HPCC has already done so for us. We need to load the FASTQC module to access the executable. 
FASTQC can be run directly from the command line for several files at a time using the following 
command line:

```bash
fastqc --outdir=output_directory --extract somefile.fastq someotherfile.fastq
```

Since there are several files to quality check, we will run a master script to submit a single 
job per animal. In this case, we will send 12 qsub texts to SLURM, and each script will generate 
four summary files, two for R1 and two for R2 fastq files (one for each sequencing run). 
First, let us create the master script:

```bash
nano $HOME/RNAseq_Pipeline/FastQC_raw.sh
```

> Copy the following script and paste in the terminal editor window.


```bash
#==============================================================================
#   File: FastQC_raw.sh
#   Directory code: $HOME/RNAseq_Pipeline/Quality/Raw
#   Date: January 6, 2020
#   Description: Generate FastQC reports for raw sequence reads.
#   Run: bash FastQC_raw.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/RAW/*fastq
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Quality/Raw
#==============================================================================

### Directory input files
raw=$SCRATCH/RAW

### Directory output files
## Create output directory
mkdir $HOME/RNAseq_Pipeline/Quality
mkdir $HOME/RNAseq_Pipeline/Quality/Raw

# Output directory
out=$HOME/RNAseq_Pipeline/Quality/Raw

## Create qstat directory for SLURM input/output files
mkdir $HOME/RNAseq_Pipeline/Quality/Raw/qstat

# SLURM input/output directory
qstat=$HOME/RNAseq_Pipeline/Quality/Raw/qstat

# Move bash script to Raw directory
mv $HOME/RNAseq_Pipeline/FastQC_raw.sh $out

# Animal IDs
cd $raw
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

# Change to raw fastq directory
cd '$raw'

# Generate quality reports
fastqc --outdir='$out' '${anim[$i]}'*

# qstat
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}_FastQC.qsub

# Submit job to hpcc
cd $qstat
sbatch ${anim[$i]}_FastQC.qsub

done
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.

Run the FASTQC master script with `bash`:

```bash
bash $HOME/RNAseq_Pipeline/FastQC_raw.sh
```
The FASTQC master will create three directories in your home space. We will recheck 
quality scores again after quality filtering. The `Quality` directory will house all 
quality check files. The `Raw` directory will contain only files of the raw data quality 
assessment that we just submitted. The third directory, `qstat` is for the job scripts 
generated with the master FASTQ executable. Let us take a look at one of the submitted 
scripts:

```bash
cd $HOME/RNAseq_Pipeline/Quality/Raw/qstat
cat `ls *.qsub | head -1`
```

The asterisk `*` metacharacter in Unix is used to match any character. In this case, we are using it to point to all 
the qsub scripts and selecting to read only the first one.

Check the status of the submitted jobs by using the show jobs shortcut:

```bash
sq
```

For now, we wait for the jobs to finish. Only when the `sq` displays zero running jobs 
can we proceed to the next step.


### Check FASTQC output with MultiQC

Before checking the quality summaries for the raw sequence, files let us first make sure all 
our jobs ran without errors or warnings. Change to the `qstat` directory and pick one of the 
output files to read using `cat`.

```bash
cd $HOME/RNAseq_Pipeline/Quality/Raw/qstat
cat `ls *FastQC.o* | head -1`
```

To review all the output files, check for errors, and obtain run summary, use my in-house 
check jobs script. Within the qstat directory run the following command:

```bash
bash $HOME/.checkJobs
```

This script will look for errors or warnings in job outputs and print it to the command line. 
A job summary is added to each output file that includes used memory and compute time. 
This information is useful to inform future resource usage. Look over these summaries by 
reviewing the same output file as before using the `cat` command.

```bash
cat `ls *FastQC.o* | head -1`
```

Running FASTQC would generate several HTML summary files, one per fastq file. To look over the 
HTML files, you need first to download the files to your computer desktop. Open a new terminal window 
and change to your Desktop directory (If you are using a MAC or have a terminal installed). 
Copy the following command, replacing `username` with your netid, to download the FASTQ files to the desktop. 
Notice we are using `-pr` when calling `scp` to download files recursively. 

You can also transfer files using the GUI [FileZilla](https://wiki.hpcc.msu.edu/display/ITH/File+transfer).

```bash
# This command will copy the ~/RNAseq_Pipeline/Quality/Raw directory to your computer desktop:
cd Desktop
scp -pr username@gateway.hpcc.msu.edu:/mnt/username/RNAseq_Pipeline/Quality/Raw .
```

Enter your password when prompted. After the files have downloaded, move to the Raw directory 
and review a few of the HTML file. Refer to [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
for a description of these files and comparisons of GOOD versus BAD quality summaries. 
We do not want to waste time opening every file to asses quality. A better approach is to merge 
all the FASTQC results into a single file that we can compare across samples. To do so, we use 
the [MultiQC](https://multiqc.info/) software that concatenates FASTQC results in a single file. 
MultiQC can also create summaries for alignment statistics, but for now, we will focus on quality 
scores. HPCC does not have this software, but we can install it using a virtual environment. 
We do not need to run this software per file; instead, we point the executable to the directory 
containing the data of interest. The following script runs MultiQC to concatenate all FASTQC files.

```bash
nano $HOME/RNAseq_Pipeline/multiQC_raw.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH -J 'MultiQC'
#SBATCH -o 'MultiQC.o%j'

#==============================================================================
#   File: multiQC_raw.sh
#   Directory code: $HOME/RNAseq_Pipeline/Quality/Raw/MultiQC
#   Date: January 6, 2020
#   Description: Concatenate FASTQC quality reports into a single file.
#   Run: sbatch multiQC_raw.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $HOME/RNAseq_Pipeline/Quality/Raw
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Quality/Raw/MultiQC
#==============================================================================

### Directory input files
Qrep=$HOME/RNAseq_Pipeline/Quality/Raw

### Directory output files
## Create MultiQC output directory
mkdir $HOME/RNAseq_Pipeline/Quality/Raw/MultiQC

# Output directory
QC=$HOME/RNAseq_Pipeline/Quality/Raw/MultiQC

#' Create virtual environment
cd $QC
virtualenv QC
source QC/bin/activate

#' Install MultiQC in virtual environment
pip install multiqc

# Generate a summary file
multiqc *

# Move bash script to work directory
mv $HOME/RNAseq_Pipeline/multiQC_raw.sh $QC

# Job details
scontrol show job $SLURM_JOB_ID
```

> Save the script by pressing `Ctrl o` and `enter`. To exit press `Ctrl x`.


To run MultiQC, copy the above script to your `Raw` directory and then submit the script to 
the scheduler with `sbatch`.

```bash
sbatch $HOME/RNAseq_pipeline/multiQC_raw.sh
```

MultiQC will generate an HTML file with the concatenated FASTQC reportes.
Download this HTML file to your computer.

For more information on the MultiQC single report, summaries see [MultiQC web page](https://multiqc.info/). 
A few things you can review in the aggregated report include sequence quality histograms, sequence length 
distribution, and overrepresented sequences. Look for any abnormalities in the data. Do not be too concerned 
with adapter contamination, unless it is excessive, and a few low-quality basses, we will process these artifacts 
out in the following steps. 


I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.

