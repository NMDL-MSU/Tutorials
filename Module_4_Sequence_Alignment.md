---
title: RNA Sequence Alignment
author: Deborah Velez-Irizarry
date: Updated Jan 3 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---



### Description
In this tutorial, we will cover how to download a reference genome, 
build a reference index with Bowtie2 and short read alignment with HISAT.

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
> Assembly report  
> Assembly statistics  
> Assembled chromosome FASTA files  

**EquCab3**

The EquCab3.sh bash script will obtain the following horse reference genome files from NCBI:

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/

- README_Equus_caballus_annotation_release_103
- GCF_002863925.1_EquCab3.0_assembly_report.txt 
- GCF_002863925.1_EquCab3.0_assembly_stats.txt  
- GCF_002863925.1_EquCab3.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/

```bash
cd ~/RNAseq_Pipeline
nano EquCab3.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#######################################################################
# Script: EquCab3.sh
# Directory Script: /mnt/home/netid/RNAseq_Pipeline/Reference
# Date: January 3, 2020
#
# Description: Obtain reference genome from ncbi and concatenate 
#    chromosomal fasta file to a single refence file
#
# Input:
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#   Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/
#   GCF_002863925.1_EquCab3.0_assembly_structure/Primary_Assembly/
#   assembled_chromosomes/FASTA/
#
# Output:
#   /mnt/home/netid/RNAseq_Pipeline/Reference/chr*.fna.gz
#   /mnt/home/netid/RNAseq_Pipeline/Reference/Equus_caballus_3.0.fna.gz
#
# Run: bash EquCab3.sh
######################################################################

# Directory for horse reference genome
mkdir /mnt/home/netid/RNAseq_Pipeline/Reference
dir=/mnt/home/netid/RNAseq_Pipeline/Reference
cd $dir

# Move bash script to work directory
mv ~/RNAseq_Pipeline/EquCab3.sh $dir

# Download Annotation Reports
ftp=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/
wget -r -np -nH --cut-dirs=7 $ftp/README_Equus_caballus_annotation_release_103
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_assembly_report.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002863925.1_EquCab3.0_assembly_stats.txt

# Downlad assembled chromosomes
wget -r -np -nH --cut-dirs=11 $ftp/GCF_002863925.1_EquCab3.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/*

# Unzip Files
gunzip $dir/*fna.gz

# Create list of ordered fasta files to concatenate
# Order Chromosomes
ls | cut -f1 -d. | sed 's/chr//' | sort -n > $dir/chr.txt
tail -31 chr.txt > tmp; mv tmp chr.txt
echo X >> chr.txt
chr=(`cat $dir/chr.txt`)

# List of fasta files to concatenate
for ((i=0; i<${#chr[@]} ; i++ )) do
    echo chr${chr[$i]}.fna >> $dir/equcab.txt
done

# Complete Horse genome 3.0 build
cat $(grep -v '^#' $dir/equcab.txt) > $dir/Equus_caballus_3.0.fna

# Remove temporary files
rm chr.txt equcab.txt

# Zip all fasta files
gzip *.fna
```

> Replace `netid` with your username and run job interactively.

```bash
sed -i 's/netid/username/g' EquCab3.sh
bash EquCab3.sh
```

**BosTau9**

The BosTau9.sh bash script will obtain the following cattle reference genome files from NCBI:

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/

- README_Bos_taurus_annotation_release_106  
- GCF_002263795.1_ARS-UCD1.2_assembly_report.txt                              
- GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt
- GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/

```bash
cd ~/RNAseq_Pipeline
nano BosTau9.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#######################################################################
# Script: BosTau9.sh
# Directory Script: /mnt/home/netid/RNAseq_Pipeline/Reference
# Date: January 3, 2020
#
# Description: Obtain reference genome from ncbi and concatenate 
#    chromosomal fasta file to a single refence file
#
# Input:
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#   Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/
#   GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/
#   assembled_chromosomes/FASTA/
#
# Output:
#   /mnt/home/netid/RNAseq_Pipeline/Reference/chr*.fna.gz
#   /mnt/home/netid/RNAseq_Pipeline/Reference/Bos_taurus_9.fna.gz
#
# Run: bash BosTau9.sh
######################################################################

# Directory for horse reference genome
mkdir /mnt/home/netid/RNAseq_Pipeline/Reference
dir=/mnt/home/netid/RNAseq_Pipeline/Reference
cd $dir

# Move bash script to work directory
mv ~/RNAseq_Pipeline/bosTau9.sh $dir

# Download Annotation Reports
ftp=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/
wget -r -np -nH --cut-dirs=7 $ftp/README_Bos_taurus_annotation_release_106
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt
wget -r -np -nH --cut-dirs=7 $ftp/GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt

# Downlad assembled chromosomes
wget -r -np -nH --cut-dirs=11 $ftp/GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/*

# Unzip Files
gunzip $dir/*fna.gz

# Create list of ordered fasta files to concatenate
# Order Chromosomes
ls | cut -f1 -d. | sed 's/chr//' | sort -n > $dir/chr.txt
tail -31 chr.txt > tmp; mv tmp chr.txt
echo X >> chr.txt
chr=(`cat $dir/chr.txt`)

# List of fasta files to concatenate
for ((i=0; i<${#chr[@]} ; i++ )) do
    echo chr${chr[$i]}.fna >> $dir/BosTau9.txt
done

# Complete cattle genome 9 build
cat $(grep -v '^#' $dir/BosTau9.txt) > $dir/Bos_taurus_9.fna

# Remove temporary files
rm chr.txt BosTau9.txt

# Zip all fasta files
gzip *.fna
```

> Replace `netid` with your username and run job interactively.

```bash
sed -i 's/netid/username/g' BosTau9.txt.sh
bash BosTau9.txt.sh
```

### Index reference genome 

We need to index the reference genome to efficiently map the short RNA-seq reads. 
This will be accomplished with Bowtie2, an ultrafast and memory-efficient tool for 
aligning sequencing reads to long reference sequences. For more information on
Bowtie2 please review the following references and webcite.

**Bowtie2 References**

Langmead B, Wilks C, Antonescu V, Charles R. [Scaling read aligners to hundreds of threads on general-purpose processors](https://www.ncbi.nlm.nih.gov/pubmed/30020410). Bioinformatics. 2018 Jul 18. doi: 10.1093/bioinformatics/bty648.
Langmead B, Salzberg SL. [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923). Nature Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.

[Bowtie2 Webcite](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)

```bash
cd ~/RNAseq_Pipeline
nano build_reference_index_bowtie2.sh
```

> Copy the following script and paste in the terminal editor window.

```bash
#==============================================================================
#   File: build_reference_index_bowtie2.sh
#   Directory code: /mnt/home/netid/RNAseq_Pipeline/Reference/bowtie2
#   Date: January 3, 2020
#   Description: Index reference genome for use in the alignment of RNA-seq
#         reads.
#   Run: bash build_reference_index_bowtie2.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       /mnt/home/netid/RNAseq_Pipeline/Reference
#
#   Output files to directory:
#       /mnt/home/netid/RNAseq_Pipeline/Reference/bowtie2
#==============================================================================

### Directory where reference genome is located
ref=/mnt/home/netid/RNAseq_Pipeline/Reference

### Change chromosome names in reference fna file
# Build index for chromosomes
cd $ref
grep -n Equus Equus_caballus_3.0.fna > check_chromosome_before.txt
idx=(`grep Equus Equus_caballus_3.0.fna  | cut -f1 -d' ' | sed 's/>//g'`)
chr=(`grep Equus Equus_caballus_3.0.fna | cut -f9 -d' ' | sed 's/,//'`)

# Loop through chromosomes and change name
scp  Equus_caballus_3.0.fna Equus_caballus_3.0_chr.fna
for ((i=0; i<${#idx[@]} ; i++ )) do
    sed -e 's/'${idx[$i]}'/chr'${chr[$i]}'/' Equus_caballus_3.0_chr.fna > tmp; mv tmp Equus_caballus_3.0_chr.fna
done
grep -n Equus Equus_caballus_3.0_chr.fna > check_chromosome_after.txt

# **Review the check files to make sure the chromosome names are correct**

### Create bowtie2 Work directory
mkdir $ref/bowtie2
index=$ref/bowtie2

# Move reference genome to bowtie directory
mv Equus_caballus_3.0_chr.fna $index/Index_EquCab3.fa
mv check_chromosome_*.txt $index

### Save index batch script and submit to hpcc
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=25G
#SBATCH -J EquCab3
#SBATCH -o '$index'/EquCab3.o%j

module load ncurses
module load ifort/2017.4.196-GCC-6.4.0-2.28
module load impi
module load Bowtie2/2.3.4

# Work directory
index='$index'

# Index Genome
cd $index
bowtie2-build Index_EquCab3.fa Index_EquCab3

# qstat
scontrol show job $SLURM_JOB_ID' > $index/index_EquCab3.qsub

# Move bash script to work directory
mv ~/RNAseq_Pipeline/build_reference_index_bowtie2.sh $index

# Submit script to hpcc
cd $index
sbatch index_EquCab3.qsub
```

> Replace `netid` with your username and run job interactively.

```bash
sed -i 's/netid/username/g' build_EquCab3_index_bowtie2.sh
bash build_EquCab3_index_bowtie2.sh
```

Notice that before we index the reference genome we update the chromosome names.
This step is added to avoid the use of chromosome identifiers in the subsequent steps.
Before procceding to the next step check that the bowtie2 job finished with no errors.

```bash
cd ~/RNAseq_Pipeline/Reference/bowtie2
bash ~/RNAseq_Pipeline/.checkJobs
```

### Alignment


