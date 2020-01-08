
---
title: RNA Sequence Alignment
author: Deborah Velez-Irizarry
date: Updated Jan 8 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---

```{r echo = FALSE, message = FALSE}
opts_chunk$set(tidy=TRUE)
```

### Description
In this tutorial, we will cover how to download a reference genome, 
build a reference index with HISAT2 build and short read alignment with HISAT2.

![](https://user-images.githubusercontent.com/44003875/71802562-00e62b00-302c-11ea-9584-bef4ec248f61.png)


### Connect to HPC system
If you are using Ubuntu, open up the terminal with `Ctrl`+`Alt`+`T`. 
From the command line, log in to HPC using ssh. 

```bash
ssh -YX username@gateway.hpcc.msu.edu
```

If you are using the remote desktop environment login to your terminal 
through the Web-based remote desktop: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)

### Obtain Referece Genome
We will use NCBI to download your reference genome. Below are two examples, 
the first will download the most horse genome (EquCab3) and the second the
cattle genome (ARS-UCD1.2). If you wish to download a different reference 
go to [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/)
and note the path to your reference of interest.

The files we will download from NCBI in this tutorial are the following:  
> Annotation release information. 
> Assembly report and statistics   
> Reference FASTA and GTF (annotation) files  

**EquCab3**

The EquCab3.sh bash script will obtain the following horse reference genome files from NCBI:

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/

- README_Equus_caballus_annotation_release_103
- GCF_002863925.1_EquCab3.0_assembly_report.txt 
- GCF_002863925.1_EquCab3.0_assembly_stats.txt  
- GCF_002863925.1_EquCab3.0_genomic.fna.gz
- GCF_002863925.1_EquCab3.0_genomic.gtf.gz

```bash
nano $HOME/RNAseq_Pipeline/EquCab3.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#######################################################################
# Script: EquCab3.sh
# Directory Script: $HOME/RNAseq_Pipeline/Reference
# Date: January 6, 2020
#
# Description: Obtain reference genome from ncbi and concatenate
#    chromosomal fasta file to a single refence file
#
# Input:
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#   Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0
#
# Output Directory:
#   $HOME/RNAseq_Pipeline/Reference
#
# Output Files:
#	README_Equus_caballus_annotation_release_103
#	GCF_002863925.1_EquCab3.0_assembly_report.txt
#	GCF_002863925.1_EquCab3.0_assembly_stats.txt
#	GCF_002863925.1_EquCab3.0_genomic.fna.gz 
#	GCF_002863925.1_EquCab3.0_genomic.gtf.gz
#	GenBankAccn.txt
#	EquCab3.fa
#	EquCab3.gtf
#
# Run: batch EquCab3.sh
######################################################################

# Directory for horse reference genome
mkdir $HOME/RNAseq_Pipeline/Reference
dir=$HOME/RNAseq_Pipeline/Reference
cd $dir

### Write qsub script and submit to HPCC
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH -J EquCab3
#SBATCH -o EquCab3.o%j

### Download annotation reports and reference genome
cd '$dir'
ftp=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0
wget -r -np -nH --cut-dirs=7 $ftp/README_Equus_caballus_annotation_release_103
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_assembly_report.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_assembly_stats.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_genomic.fna.gz
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_genomic.gtf.gz   

# Unzip Files
gunzip *.gz

# Sequence name gene bank accession ID
cnt=(`cat *assembly_report.txt | wc -l`)
hd=(`grep -n '"'"'#'"'"' *assembly_report.txt | tail -1 | cut -f1 -d:`)
tail -$(($cnt - $hd)) *assembly_report.txt | cut -f1,4,7 > GenBankAccn.txt

# Create copy
scp -p *genomic.fna EquCab3.fa
scp -p *genomic.gtf EquCab3.gtf

# Zip Files
gzip *genomic*

### Change chromosome names in fa/gtf file
nm=(`ls *.fa | cut -f1 -d.`)

# Loop through chromosomes and change name
idx=(`grep Chromosome GenBankAccn.txt | cut -f3`)
chr=(`grep Chromosome GenBankAccn.txt | cut -f1`)
for ((i=0; i<${#idx[@]} ; i++ )) do
    sed -e '"'"'s/'"'"'${idx[$i]}'"'"'/'"'"'${chr[$i]}'"'"'/'"'"' $nm.fa > tmp; mv tmp $nm.fa
    sed -e '"'"'s/'"'"'${idx[$i]}'"'"'/'"'"'${chr[$i]}'"'"'/'"'"' $nm.gtf > tmp; mv tmp $nm.gtf
done' > $dir/EquCab3.qsub

# Submit script to hpcc
cd $dir
sbatch EquCab3.qsub

### Move bash script to work directory
mv $HOME/RNAseq_Pipeline/EquCab3.sh $dir
```

> Run job interactively.

```bash
bash $HOME/RNAseq_Pipeline/EquCab3.sh
```

**BosTau9**

The BosTau9.sh bash script will obtain the following cattle reference genome files from NCBI:

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/

- README_Bos_taurus_annotation_release_106  
- GCF_002263795.1_ARS-UCD1.2_assembly_report.txt                              
- GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt
- GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
- GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz

```bash
nano $HOME/RNAseq_Pipeline/BosTau9.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#######################################################################
# Script: BosTau9.sh
# Directory Script: $HOME/RNAseq_Pipeline/Reference
# Date: January 6, 2020
#
# Description: Obtain reference genome from ncbi and concatenate
#    chromosomal fasta file to a single refence file
#
# Input:
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#   Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/
#
# Output Directory:
#   $HOME/RNAseq_Pipeline/Reference
#
# Output Files:
#   README_Bos_taurus_annotation_release_106
#	GCF_002263795.1_ARS-UCD1.2_assembly_report.txt
#	GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt
#	GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz 
#	GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz
#	BosTau9.fa
#	GenBankAccn.txt
#
# Run: batch BosTau9.sh
######################################################################

# Directory for horse reference genome
mkdir $HOME/RNAseq_Pipeline/Reference
dir=$HOME/RNAseq_Pipeline/Reference
cd $dir

### Write qsub script and submit to HPCC
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH -J BosTau9
#SBATCH -o BosTau9.o%j

### Download annotation reports and reference genome
ftp=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/
wget -r -np -nH --cut-dirs=7 $ftp/README_Bos_taurus_annotation_release_106
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz  
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz

# Unzip Files
gunzip *.gz

# Sequence name gene bank accession ID
cnt=(`cat *assembly_report.txt | wc -l`)
hd=(`grep -n '"'"'#'"'"' *assembly_report.txt | tail -1 | cut -f1 -d:`)
tail -$(($cnt - $hd)) *assembly_report.txt | cut -f1,4,7 > GenBankAccn.txt

# Create copy
scp -p *genomic.fna BosTau9.fa
scp -p *genomic.gtf BosTau9.gtf

# Zip Files
gzip *genomic*

### Change chromosome names in fa/gtf file
nm=(`ls *.fa | cut -f1 -d.`)

# Loop through chromosomes and change name
idx=(`grep Chromosome GenBankAccn.txt | cut -f3`)
chr=(`grep Chromosome GenBankAccn.txt | cut -f1`)
for ((i=0; i<${#idx[@]} ; i++ )) do
    sed -e '"'"'s/'"'"'${idx[$i]}'"'"'/chr'"'"'${chr[$i]}'"'"'/'"'"' $nm.fa > tmp; mv tmp $nm.fa
    sed -e '"'"'s/'"'"'${idx[$i]}'"'"'/chr'"'"'${chr[$i]}'"'"'/'"'"' $nm.gtf > tmp; mv tmp $nm.gtf
done' > $dir/BosTau9.qsub

# Submit script to hpcc
cd $dir
sbatch BosTau9.qsub

### Move bash script to work directory
mv $HOME/RNAseq_Pipeline/BosTau9.sh $dir
```

> Run job interactively.

```bash
bash $HOME/RNAseq_Pipeline/BosTau9.txt.sh
```

Notice that before we index the reference genome we update the chromosome names.
This step is added to avoid the use of chromosome identifiers in the subsequent steps.


### Index reference genome 

We need to index the reference genome to efficiently map the short RNA-seq reads. 
This will be accomplished with HISAT2 build, an ultrafast and memory-efficient tool for 
aligning sequencing reads to long reference sequences. For more information on
HISAT2 please review the following references and webcite.

**HISAT2 References**

Kim, Daehwan and Langmead, Ben and Salzberg, Steven L (2015). 
[HISAT: a fast spliced aligner with low memory requirements.](https://www.ncbi.nlm.nih.gov/pubmed/25751142) 
INature Methods, 12(4):357-360.

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL (2016).
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown.](https://www.ncbi.nlm.nih.gov/pubmed/27560171)
Nature Protocols, 11(9):1650-1667

[HISAT2 Manual](https://ccb.jhu.edu/software/hisat2/manual.shtml)

```bash
nano $HOME/RNAseq_Pipeline/build_reference_index_hisat2.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: build_reference_index_hisat2.sh
#   Directory code: $HOME/RNAseq_Pipeline/Reference/Index
#   Date: January 6, 2020
#   Description: Index reference genome for use in the alignment of RNA-seq
#         reads.
#   Run: bash build_reference_index_hisat2.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $HOME/RNAseq_Pipeline/Reference
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/Reference/Index
#==============================================================================

### Directory where reference genome is located
ref=$HOME/RNAseq_Pipeline/Reference

### Build index for chromosomes

#Create Index directory
mkdir $ref/Index
index=$ref/Index

# Reference genome file
cd $ref
nm=(`ls *.fa | cut -f1 -d.`)

### Write qsub script and submit to HPCC
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=04:00:00
#SBATCH --mem=25G
#SBATCH -J Index
#SBATCH -o '$index'/Index.o%j

# Load required modules
module load SAMtools/1.9
module load hisat2/2.1.0

# Module List
module list

# Index Genome
cd '$index'
hisat2-build -p 10 '$ref/$nm'.fa Index_'$nm'

# qstat
scontrol show job $SLURM_JOB_ID' > $index/Index.qsub

# Move bash script to work directory
mv $HOME/RNAseq_Pipeline/build_reference_index_hisat2.sh $index

# Submit script to hpcc
cd $index
sbatch Index.qsub
```

> Run master script.

```bash
bash $HOME/RNAseq_Pipeline/build_reference_index_hisat2.sh
```

Before procceding to the next step check that the job finished 
with no errors.

```bash
cd $HOME/RNAseq_Pipeline/Reference/Index
checkJobs
```


### Alignment of reads to reference

HISAT2 (hierarchical indexing for spliced alignment of transcripts) is an 
alignment program for mapping next-generation sequencing reads. This method 
have been described in [Nature Protocols](https://www.ncbi.nlm.nih.gov/pubmed/27560171).
We will run the alignment per animal using a master script. To ensure the jobs
finish in the time allocated we will use eight cores.

There are several options available for HISAT2, go to the [HISAT2 mannual](https://ccb.jhu.edu/software/hisat2/manual.shtml) to view
all available options.

In this tutorial we are aligning Illimina paired end reads, from fastq files, 
that have strand specific information. The options selected are:

> `-q` refering to fastq files  
> `--phred33` input qualities are ASCII chars equal to phred33 quality scores (Illumina)  
> `--rna-strandness RF` specifies the strand-specific information as `R` first strand (or reverse complement)
> and `F` the second strand (or transcript)   
> `--met-stderr` saves the metrics report to the "standard error" output file  
> `--dta-cufflinks` report strand inforamation for every spliced alignment (XS:A[+-])  
> `-p 8` parallel search threads, in this case 8  
> `-x` basename of the index for the reference genome  
> `-1` mate 1 fastq file  
> `-2` mate 2 fastq file  
> `-U` single mates (unpaired reads) fastq file  
> `-S` name for output sam file  

```bash
nano $HOME/RNAseq_Pipeline/hisat2_align_reads.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: hisat2_align_reads.sh
#   Directory code: $HOME/RNAseq_Pipeline/HISAT2
#   Date: January 6, 2020
#   Description: Run the splice junction mapper using the stranded
#                protocol with HISAT2 on the adapter and quality trimmed data.
#                A folder will be created per sample.
#   Run: bash hisat2_align_reads.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/Trimmomatic/PE
#       $SCRATCH/Trimmomatic/SE
#       /$HOME/RNAseq_Pipeline/Reference/Index
#
#   Output files to directory:
#       $SCRATCH/HISAT2
#==============================================================================

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Work Directory
mkdir $HOME/RNAseq_Pipeline/HISAT2
dir=$HOME/RNAseq_Pipeline/HISAT2

# Input Directory
pe=$SCRATCH/Trimmomatic/PE
se=$SCRATCH/Trimmomatic/SE
index=$HOME/RNAseq_Pipeline/Reference/Index
nm=(`ls $index/*.ht2 | cut -f1 -d. | uniq`)

# Output Directory
mkdir $SCRATCH/HISAT2
out=$SCRATCH/HISAT2

# Qstat directory
mkdir $dir/qstat
qstat=$dir/qstat

# Move script to directory
mv $HOME/RNAseq_Pipeline/hisat2_align_reads.sh $dir

# Write bash script for each animal
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=300G
#SBATCH -J '${anim[$i]}_hisat2'
#SBATCH -o '${anim[$i]}_hisat2.o%j'

# Work Directory
cd '$out'

# Load required modules
module load SAMtools/1.9
module load hisat2/2.1.0

# Module List
module list

# Merge all single reads to one file
cat '$se'/'${anim[$i]}'_R1_merged_SE.fastq > '$se'/'${anim[$i]}'_merged_SE.fastq
cat '$se'/'${anim[$i]}'_R2_merged_SE.fastq >> '$se'/'${anim[$i]}'_merged_SE.fastq

# Align reads
hisat2 -q --phred33 --rna-strandness RF --met-stderr --dta-cufflinks -p 8  -x '$nm' -1 '$pe'/'${anim[$i]}'_R1_merged_PE.fastq \
 -2 '$pe'/'${anim[$i]}'_R2_merged_PE.fastq -U '$se'/'${anim[$i]}'_merged_SE.fastq -S '${anim[$i]}'.sam

# Convert SAM to BAM format
samtools view -bS '${anim[$i]}'.sam > '${anim[$i]}'.bam

# Run statistics
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}.qsub

# Submit script to hpcc
cd $qstat
sbatch ${anim[$i]}.qsub

done
```
> Run master script.

```bash
bash $HOME/RNAseq_Pipeline/hisat2_align_reads.sh
```

Check submitted jobs `sq`. When all jobs have completed, check that the 
hisat2 jobs finished with no errors before procceding to the next step.

```bash
cd $HOME/RNAseq_Pipeline/HISAT2/qstat
checkJobs
```

### Extract alignment statistics

The alignment statistics per animal were saved to the SLURM output file.
Let us look at one of these summaries.

```bash
cd $HOME/RNAseq_Pipeline/HISAT2/qstat
cat `ls | grep -v qsub | head -1`
```

We will extract the most important statistics from this output. 
First create the master script.

```bash
nano $HOME/RNAseq_Pipeline/check_hisat2_jobs.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: check_hisat2_jobs.sh
#   Directory code: $HOME/RNAseq_Pipeline/HISAT2
#   Date: January 6, 2020
#   Description: Check that all scripts ran, produced output and
#                finished without errors.
#                Generate alignment summaries for all animals.
#                Generate bam files containing only uniquely mapped reads.
#   Run: bash check_hisat2_jobs.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/HISAT2
#
#   Output files to directory:
#       $SCRATCH/HISAT2
#
#   Output files:
#       *_stats.txt (statistics of mapper)
#       *_uniq.bam
#==============================================================================

# Input Directory
dir=$SCRATCH/HISAT2

# Output Directory
mkdir $HOME/RNAseq_Pipeline/HISAT2/qstat/summary
qstat=$HOME/RNAseq_Pipeline/HISAT2/qstat/summary

# Home directory
home=$HOME/RNAseq_Pipeline/HISAT2

# Move bash script to HISAT2 directory
mv $HOME/RNAseq_Pipeline/check_hisat2_jobs.sh $home

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Check that you have output for each file
for ((i=0; i<${#anim[@]} ; i++ )) do
ls $dir/${anim[$i]}.bam > ../temp
done; rm ../temp

# Summary statistics from alignment
for ((i=0; i<${#anim[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=25G
#SBATCH -J '${anim[$i]}_summary'
#SBATCH -o '${anim[$i]}_summary.o%j'

# Load required modules
module load SAMtools/1.9
module list

### Create summary statistics of alignment using the align_summary output from TopHat
cd '$home/qstat'
# Input
echo total_input_reads `grep reads '${anim[$i]}'_hisat2.o* | cut -f1 -d'"'"' '"'"'` > '$dir/${anim[$i]}'_stats.txt
echo paired_input_reads `grep '"'"'were paired'"'"' '${anim[$i]}'_hisat2.o* | cut -f3 -d'"'"' '"'"'` >> '$dir/${anim[$i]}'_stats.txt
echo unpaired_input_reads `grep '"'"'were unpaired'"'"' '${anim[$i]}'_hisat2.o* | cut -f3 -d'"'"' '"'"'` >> '$dir/${anim[$i]}'_stats.txt

# Aligned
echo paired_aligned_reads `grep --context=3 '"'"'were paired'"'"' '${anim[$i]}'_hisat2.o* |  tail -2 | cut  -f5 -d'"'"' '"'"' | paste -sd+ - | bc` >> '$dir/${anim[$i]}'_stats.txt
echo unpaired_aligned_reads `grep --context=3 '"'"'were unpaired'"'"' '${anim[$i]}'_hisat2.o* | tail -2 | cut  -f5 -d'"'"' '"'"' | paste -sd+ - | bc` >> '$dir/${anim[$i]}'_stats.txt

# Dropped
grep '"'"'0 times'"'"' '${anim[$i]}'_hisat2.o*| head -1 | cut -f5 -d'"'"' '"'"' > '${anim[$i]}'_drop
grep '"'"'0 times'"'"' '${anim[$i]}'_hisat2.o* | tail -1 | cut -f5 -d'"'"' '"'"' >> '${anim[$i]}'_drop
echo unaligned_reads `cat '${anim[$i]}'_drop | paste -sd+ - | bc` >> '$dir/${anim[$i]}'_stats.txt
rm drop

# Total Aligned
map=(`grep aligned '$dir/${anim[$i]}'_stats.txt | head -2 | cut -f2 -d'"'"' '"'"' | paste -sd+ - | bc`)
tot=(`cat '$dir/${anim[$i]}'_stats.txt | cut -f2 -d'"'"' '"'"' | head -1`)
echo total_aligned_reads `echo $map` >> $dir/'${anim[$i]}'_stats.txt
echo alignment_rate '"'"'0'"'"'`echo "scale=6 ; ($map / $tot)" | bc` >> '$dir/${anim[$i]}'_stats.txt

# Total Uniquely Aligned
grep '"'"'exactly 1 time'"'"' '${anim[$i]}'_hisat2.o* | cut -f5 -d'"'"' '"'"' | head -1 > '${anim[$i]}'_uni
grep '"'"'exactly 1 time'"'"' '${anim[$i]}'_hisat2.o* | cut -f5 -d'"'"' '"'"' | tail -1 >> '${anim[$i]}'_uni
echo total_uniquely_aligned_reads `cat '${anim[$i]}'_uni | paste -sd+ - | bc` >> '$dir/${anim[$i]}'_stats.txt
map=(`cat '${anim[$i]}'_uni | paste -sd+ - | bc`)
echo unique_alignment_rate '"'"'0'"'"'`echo "scale=6 ; ($map / $tot)" | bc` >> '$dir/${anim[$i]}'_stats.txt
rm uni

### Retain only uniquely aligned reads
(samtools view -H '$dir/${anim[$i]}'.bam; samtools view -F 2308 '$dir/${anim[$i]}'.bam | grep -w '"'"'NH:i:1'"'"') | samtools view -bS - > '$dir/${anim[$i]}'_uniq.bam

### Job details
echo 'Job Details'
scontrol show job $SLURM_JOB_ID' > $qstat/${anim[$i]}_summary.qsub

# Submit script to hpcc
cd $qstat
sbatch ${anim[$i]}_summary.qsub

done
```
> Run master script.

```bash
bash $HOME/RNAseq_Pipeline/check_hisat2_jobs.sh
```

Check submitted jobs `sq`. When all jobs have completed, check if 
jobs finished with no errors before procceding to the next step.

```bash
cd $HOME/RNAseq_Pipeline/HISAT2/qstat/summary
checkJobs
```


### HISAT2 summary file

Now that we have alignment statistics per animal we can concatinate then
to matrix form. This summary file will be saved as `summary_alignment.txt`.

```bash
nano $HOME/RNAseq_Pipeline/summary_hisat2.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: summary_hisat2.sh
#   Directory code: $HOME/RNAseq_Pipeline/HISAT2
#   Date: January 6, 2020
#   Description: Merge alignment summaries for each animal to a single file.
#   Run: bash summary_hisat2.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       $SCRATCH/HISAT2
#
#   Output files to directory:
#       $HOME/RNAseq_Pipeline/HISAT2/summary_alignment.txt
#==============================================================================

# Input Directory
dir=$SCRATCH/HISAT2

# Output Directory
out=$HOME/RNAseq_Pipeline/HISAT2

# Move bash script to HISAT2 directory
mv $HOME/RNAseq_Pipeline/summary_hisat2.sh $out

# Animal IDs
cd $SCRATCH/RAW
anim=(`ls *.fastq | cut -f1 -d_ | uniq`)

# Summary statistics from alignment
cd $dir
for ((i=0; i<${#anim[@]} ; i++ )) do

    # Add column names to file
    if [ $i == 0 ]; then
    echo animal `cat ${anim[$i]}_stats.txt | cut -f1 -d' ' | tr "\n" " " ` > $out/summary_alignment.txt
    echo "" >> $out/summary_alignment.txt
    fi

    # Add alignment statisticas per animal to file
    echo ${anim[$i]} `cat ${anim[$i]}_stats.txt | cut -f2 -d' ' | tr "\n" " " ` >> $out/summary_alignment.txt
done
```

> Run script interactively.

```bash
bash $HOME/RNAseq_Pipeline/summary_hisat2.sh
```

Let us look at the summary file

```bash
head $HOME/RNAseq_Pipeline/HISAT2/summary_alignment.txt
```

I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.



