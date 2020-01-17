---
title: RNA-Seq Pipeline Statistics
author: Deborah Velez-Irizarry
date: Mon Oct 22 14:12:28 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Alignment statistics for RER Thoroughbred after mapping reads to the reference genome EquCab3  
  
***  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis/Descriptive_Analysis_RER_Thoroughbred.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;Trimmomatic/trimmomatic_rst.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/paired/rst_quality_trim_paired.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R1.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R2.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Tophat/summary_alignment.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Tophat/uniq_depth.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/summary_htseq.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/paired_counts.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/single_counts.txt  
>&nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt
  
**Output files:**  
  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis_KER_Glycogen/retained_read_stats_RER_Thoroughbred.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/htseq_counts_RER_Thoroughbred.txt 
 
***  
### Code  
Clear Environment


```r
rm(list=ls())
```

**Session Information**


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.20
## 
## loaded via a namespace (and not attached):
## [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     stringi_1.2.3  
## [5] stringr_1.3.1   evaluate_0.10.1
```

**Input Directory**


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
```

### Animal IDs


```r
animID <- rownames(read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"), 
    header=TRUE, row.names=1))[1:23]
animID
```

```
##  [1] "12401" "12402" "12403" "12610" "12611" "12612" "12613" "12614"
##  [9] "12616" "12617" "12618" "12619" "12620" "12621" "12622" "12910"
## [17] "12913" "12915" "12916" "12917" "12918" "12921" "12924"
```

### Load RNA-Seq Pipeline Statistics:  
**Trimmomatic**


```r
adpt <- read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"), 
    header=TRUE, row.names=1)
adpt <- adpt[animID,]
dim(adpt)
```

```
## [1] 23  5
```

**Condetri**  
> Paired reads


```r
cond.PE <- read.table(paste(dir,"Condetri/paired", 
    "rst_quality_trim_paired.txt", sep="/"), header=TRUE, row.names=1)
cond.PE <- cond.PE[animID,]
dim(cond.PE)
```

```
## [1] 23  6
```

> Single Reads: R1


```r
cond.SR1 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R1.txt", sep="/"), header=TRUE, row.names=1)
cond.SR1 <- cond.SR1[paste(animID, "R1", sep="_"),]
dim(cond.SR1)
```

```
## [1] 23  6
```

> Single Reads: R2


```r
cond.SR2 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R2.txt", sep="/"), header=TRUE, row.names=1)
cond.SR2 <- cond.SR2[paste(animID, "R2", sep="_"),]
dim(cond.SR2)
```

```
## [1] 23  6
```

**Tophat**


```r
tophat <- read.table(paste(dir,"Tophat", "summary_alignment.txt", sep="/"), 
    header=TRUE, row.names=1)
tophat <- tophat[animID,]
dim(tophat)
```

```
## [1] 23  7
```

**Depth**  
> Uniquely mapped reads (includes paired and single reads)


```r
uniq.depth <- read.table(paste(dir,"Tophat", "uniq_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.depth <- uniq.depth[animID,]
dim(uniq.depth)
```

```
## [1] 23  3
```

> Uniquely mapped reads per chromosome (includes paired and single reads)


```r
uniq.chr.depth <- read.table(paste(dir,"Tophat", "uniq_chr_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.chr.depth <- uniq.chr.depth[animID,]
dim(uniq.chr.depth)
```

```
## [1] 23 33
```

**HTSeq**  
> Summary HTSeq table


```r
sumHTSeq <- read.table(paste(dir, "HTSeq", "summary_htseq.txt", sep="/"), 
    header=TRUE, row.names=1)
sumHTSeq <- sumHTSeq[animID,]
dim(sumHTSeq)
```

```
## [1] 23  4
```

> Paired read counts table


```r
paired <- read.table(paste(dir, "HTSeq", "paired_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
paired <- paired[,paste("X", animID, "_paired", sep="")]
dim(paired)
```

```
## [1] 37875    23
```

> Single read count table


```r
single <- read.table(paste(dir, "HTSeq", "single_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
single <- single[,paste("X", animID, "_single", sep="")]
dim(single)
```

```
## [1] 37875    23
```

> Annotation


```r
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation", "annotation.txt", sep="/"),
    header=TRUE, row.names=8)
```

 ### Summary Function


```r
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]
```

### Summary Trimmomatic: Adapter filtering
**Number of Input Read Paires (million reads)**


```r
startPE <- adpt$input.read.pairs*2
names(startPE) <- rownames(adpt)
summSD((startPE/2)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  30.0457  41.5475  45.4416  45.6394   6.2325  48.0077  60.0321
```

**Number of Input Read (million reads)**


```r
summSD((startPE)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  60.0914  83.0950  90.8832  91.2789  12.4650  96.0154 120.0642
```

**Number of retained paired reads after filtering adapter sequences (million)**


```r
summSD(adpt$both.surviving*2/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  47.2970  66.9871  71.1830  72.1295   9.6209  76.5366  92.9136
```

Percent retained paired reads


```r
mean((adpt$both.surviving*2) / startPE)*100
```

```
## [1] 79.05251
```

**Number of forward reads surviving without its pair (million)**


```r
summSD(adpt$fwd.only.surviving/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   6.2343   8.6663   9.3497   9.3847   1.5378   9.9730  13.3603
```

**Number of reverse reads surviving without its pair (million)**


```r
summSD(adpt$rev.only.surviving/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   0.0743   0.0830   0.0862   0.0922   0.0147   0.0986   0.1261
```

**Number of single reads surviving (million)**


```r
summSD((adpt$fwd.only.surviving/1e6 + adpt$rev.only.surviving/1e6))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   6.3165   8.7612   9.4493   9.4770   1.5451  10.0755  13.4580
```

Percent single reads surviving


```r
mean((adpt$fwd.only.surviving + adpt$rev.only.surviving) / startPE)*100
```

```
## [1] 10.36583
```

**Number of retained reads start after adapter trimming**


```r
trimm.out <- rowSums(data.frame((adpt$both.surviving*2), + adpt$fwd.only.surviving + 
    adpt$rev.only.surviving))
names(trimm.out) <- rownames(adpt)
summSD(trimm.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  53.6135  74.7593  80.9556  81.6065  11.0178  86.1744 106.3716
```

**Number of dropped reads (million)**


```r
summSD(startPE - trimm.out)/1e6
```

```
##      Min.   1st Qu.    Median      Mean  Std.Dev.   3rd Qu.      Max. 
##  6.477883  8.956732  9.648329  9.672365  1.564202 10.283387 13.692619
```

Percent single reads surviving


```r
mean((startPE - trimm.out) / startPE)*100
```

```
## [1] 10.58166
```

**Percent of retained reads after adapter trimming from total number of sequenced reads**


```r
summSD(trimm.out/startPE)*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    88.30    88.96    89.52    89.42     0.61    89.77    90.56
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Trimmomatic.  
Output should be zero.


```r
sum(!rowSums(data.frame((adpt$both.surviving*2), 
    (rowSums(data.frame(adpt$fwd.only.surviving, adpt$rev.only.surviving, 
    adpt$dropped)*2)))) == startPE)
```

```
## [1] 0
```

### Summary Condetri: Quality filtering  
**Data Check:** Input reads for Condetri should be equal to output reads for Trimmomatic


```r
start.cond <- rowSums(data.frame(cond.PE$TotReads, cond.SR1$TotReads, 
    cond.SR2$TotReads))
names(start.cond) <- rownames(cond.PE)

# Do the read input for condetri reflect the read output for trimmomatic? 
# Output should be zero.
sum(!start.cond == trimm.out[rownames(start.cond)])
```

```
## [1] 0
```

```r
# Do the row names for paired and single reads condetri output match? 
# Output should be zero.
sum(!rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR1), "_"), 
    function(x) x[[1]][1])) |
    !rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR2), "_"), 
    function(x) x[[1]][1])))
```

```
## [1] 0
```

**Number of paired reads retained after quality trimming (million single)**


```r
summSD((cond.PE$PairReads)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  31.0430  41.4051  45.9820  45.9710   7.0087  48.6701  60.7334
```

Percent retained paired reads


```r
mean(cond.PE$PairReads / start.cond)*100
```

```
## [1] 56.26501
```

**Number of unpaired reads (million)**


```r
summSD((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  12.7878  19.2990  20.1552  20.2976   2.7796  21.1289  26.7602
```

Percent unpaired reads


```r
mean((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads) / start.cond)*100
```

```
## [1] 24.88845
```

**Number of dropped reads after quality filtering (million)**


```r
cond.drop <- data.frame(PE=(cond.PE$TotReads - (cond.PE$PairReads + 
    cond.PE$UnparedReads)),
    SR1=(cond.SR1$TotReads - cond.SR1$UnparedReads),
    SR2=(cond.SR2$TotReads - cond.SR2$UnparedReads))
rownames(cond.drop) <- rownames(cond.PE)
summSD(rowSums(cond.drop/1e6))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   9.7827  14.4196  15.4214  15.3378   1.9775  16.4865  18.8780
```

Percent dropped


```r
mean(rowSums(cond.drop) / start.cond)*100
```

```
## [1] 18.84653
```

**Number of retained reads after quality filtering (million)**


```r
cond.out <- rowSums(data.frame(cond.PE$PairReads, 
    cond.PE$UnparedReads, cond.SR1$UnparedReads, 
    cond.SR2$UnparedReads))
names(cond.out) <- rownames(cond.PE)
summSD(cond.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  43.8308  60.6230  66.1372  66.2687   9.3984  69.8005  87.4936
```

**Percent of retained reads after quality trimming from total number of sequenced reads**


```r
summSD(cond.out/startPE[names(cond.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    70.68    71.64    72.58    72.56     1.08    73.24    74.52
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Condetri  
Output should be zero.


```r
sum(!rowSums(data.frame(cond.out, cond.drop)) == start.cond)
```

```
## [1] 0
```

### Summary Tophat: Number of reads aligning to reference genome EquCap3  
**Data Check:** Input reads for Tophat should be equal to output reads for Condreti  


```r
# Do the read input for condetri reflect the read output for trimmomatic? Output should be zero.
start.top <- tophat$total_input_reads
names(start.top) <- rownames(tophat)
sum(sum(!start.top) == cond.out[rownames(tophat)])
```

```
## [1] 0
```

**Number of paired reads aligning to the reference**


```r
summSD(tophat$total_paired_reads/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  22.6293  30.3904  34.1375  34.0773   5.3874  36.7499  45.3542
```

Percent paired reads


```r
mean(tophat$total_paired_reads/start.top)*100
```

```
## [1] 51.35982
```

**Number single reads aligning to the reference**


```r
summSD(tophat$total_unpaired_reads/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   9.6505  14.1675  15.2521  15.3021   2.1263  16.2982  19.5170
```

Paired single


```r
mean(tophat$total_unpaired_reads/start.top)*100
```

```
## [1] 23.12915
```

**Number of unmapped reads**


```r
tophat.drop <- (start.top - (tophat$total_paired_reads + 
    tophat$total_unpaired_reads))
names(tophat.drop) <- rownames(tophat)
summSD(tophat.drop/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  11.5509  14.9961  16.1815  16.8893   3.0125  19.2630  24.9965
```

Percent unmapped


```r
mean(tophat.drop/start.top)*100
```

```
## [1] 25.51103
```

**Total number of mapped reads**


```r
tophat.out <- (tophat$total_paired_reads + tophat$total_unpaired_reads)
names(tophat.out) <- rownames(tophat)
summSD(tophat.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  32.2799  44.1013  49.3542  49.3794   7.2969  52.9849  64.0908
```

**Mapping percent from total input reads**


```r
summSD(tophat$read_mapping_rate_percent)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  70.3937  72.3278  73.8049  74.4890   2.7864  76.3087  80.4436
```

Double check you get the same percent of mapped reads


```r
mean(tophat.out/start.top)*100
```

```
## [1] 74.48897
```

**Percent of reads retained (aligned) from total number of sequenced reads**


```r
summSD(tophat.out/startPE[names(tophat.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    49.93    52.56    53.72    54.05     2.18    55.94    58.23
```

### Samtools: Retain unique reads  
**Number of uniquely alignned reads**


```r
uniq.out <- tophat$total_unique_reads
names(uniq.out) <- rownames(tophat)
summSD(uniq.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  31.3320  42.9136  47.8276  47.9838   7.0863  51.4442  62.2212
```

Percent unique reads


```r
mean(uniq.out/tophat.out)*100
```

```
## [1] 97.17546
```

**Number of dropped reads**


```r
summSD(tophat.out - uniq.out)/1e6
```

```
##      Min.   1st Qu.    Median      Mean  Std.Dev.   3rd Qu.      Max. 
## 0.9478680 1.2141915 1.4405980 1.3956428 0.2344826 1.5911035 1.8695880
```

Percent dropped


```r
mean((tophat.out - uniq.out)/tophat.out)*100
```

```
## [1] 2.824538
```

**Number of uniquely algned and properly paired reads**


```r
summSD(tophat$total_unique_properly_paired/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  19.9543  26.8049  29.9144  30.1676   4.9391  32.6979  40.6688
```

**Total number of reads for HTSeq from total number of sequenced reads**


```r
summSD(uniq.out/startPE[names(uniq.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    48.62    51.13    52.24    52.52     2.05    54.37    56.44
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Tophat.  
Output should be zero.


```r
sum(!rowSums(data.frame(tophat.out, tophat.drop)) == start.top)
```

```
## [1] 0
```

```r
### Summary Depth per animal:  
```

Ratio of total nucleotides mapped to total nucleotides sequenced per animal 


```r
depth <- uniq.depth$Average_coverage
names(depth) <- rownames(uniq.depth)
summSD(depth)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  44.4631  53.9493  58.9809  60.0202   7.7533  64.1347  74.5453
```

Summary of total nucleotides mapped per chromosome for EquCap3


```r
ref <- uniq.chr.depth["EquCab3",]
depth.chr <- uniq.chr.depth[-64,]
round(t(apply(depth.chr, 2, summSD)/1e6), 2)
```

```
##         Min. 1st Qu. Median  Mean Std.Dev. 3rd Qu.  Max.
## chr1    2.97    3.99   4.24  4.30     0.63    4.64  5.67
## chr2    0.75    1.16   1.28  1.31     0.22    1.46  1.72
## chr3    0.59    0.93   1.00  1.04     0.18    1.16  1.40
## chr4    0.62    0.92   1.08  1.05     0.17    1.13  1.35
## chr5    0.92    1.43   1.59  1.64     0.29    1.82  2.24
## chr6    1.92    2.51   2.87  2.83     0.48    3.09  3.96
## chr7    1.06    1.43   1.63  1.64     0.29    1.85  2.36
## chr8    0.67    1.02   1.11  1.13     0.18    1.24  1.46
## chr9    0.42    0.64   0.71  0.71     0.12    0.77  0.94
## chr10   2.14    2.72   3.01  2.97     0.43    3.18  3.84
## chr11   1.94    2.80   3.10  3.23     0.64    3.67  4.42
## chr12   1.49    1.86   2.13  2.13     0.37    2.39  2.88
## chr13   2.05    2.44   2.78  2.81     0.52    3.19  4.07
## chr14   1.44    2.06   2.21  2.25     0.33    2.44  3.02
## chr15   0.49    0.77   0.88  0.88     0.16    0.97  1.17
## chr16   0.55    0.88   0.94  0.98     0.17    1.09  1.29
## chr17   0.37    0.53   0.60  0.61     0.11    0.66  0.85
## chr18   0.73    1.34   1.58  1.53     0.39    1.75  2.36
## chr19   0.26    0.42   0.47  0.48     0.09    0.53  0.64
## chr20   0.38    0.62   0.70  0.70     0.12    0.79  0.93
## chr21   0.29    0.42   0.47  0.48     0.08    0.53  0.62
## chr22   0.86    1.09   1.26  1.27     0.23    1.39  1.88
## chr23   0.17    0.27   0.31  0.31     0.06    0.34  0.42
## chr24   0.32    0.49   0.53  0.56     0.10    0.61  0.76
## chr25   0.53    0.69   0.75  0.76     0.11    0.85  0.96
## chr26   0.13    0.21   0.23  0.24     0.04    0.26  0.31
## chr27   4.80    5.61   6.22  6.66     1.46    7.75 10.49
## chr28   0.76    1.00   1.11  1.12     0.18    1.25  1.56
## chr29   0.12    0.20   0.23  0.23     0.05    0.26  0.30
## chr30   0.25    0.38   0.46  0.44     0.08    0.49  0.56
## chr31   0.59    0.72   0.78  0.87     0.19    0.99  1.22
## chrX    0.48    0.73   0.83  0.84     0.14    0.91  1.12
## chrAll 31.33   42.91  47.83 47.98     7.09   51.44 62.22
```

Ratio of total nucleotides mapped to total sequenced per animal


```r
ratio.chr <- do.call(rbind, lapply(colnames(depth.chr), 
    function(x) summSD(depth.chr[,x]/depth.chr[,"chrAll"])* 100))
row.names(ratio.chr) <- colnames(depth.chr)
ratio.chr
```

```
##          Min. 1st Qu. Median   Mean Std.Dev. 3rd Qu.   Max.
## chr1     8.48    8.63   8.83   8.98     0.42    9.33   9.90
## chr2     2.34    2.67   2.75   2.72     0.15    2.83   2.99
## chr3     1.84    2.07   2.17   2.16     0.15    2.26   2.39
## chr4     1.93    2.14   2.19   2.18     0.10    2.23   2.38
## chr5     2.92    3.27   3.43   3.40     0.23    3.58   3.83
## chr6     5.04    5.64   5.98   5.89     0.36    6.15   6.44
## chr7     3.11    3.27   3.37   3.40     0.20    3.54   3.79
## chr8     2.12    2.30   2.35   2.36     0.13    2.43   2.62
## chr9     1.25    1.44   1.47   1.47     0.09    1.53   1.63
## chr10    5.49    6.03   6.17   6.20     0.36    6.34   7.12
## chr11    5.58    6.17   6.47   6.71     0.74    7.09   8.53
## chr12    3.67    4.14   4.44   4.43     0.39    4.69   5.21
## chr13    4.80    5.33   5.89   5.85     0.57    6.35   6.82
## chr14    4.17    4.38   4.72   4.71     0.38    5.00   5.45
## chr15    1.52    1.74   1.85   1.83     0.13    1.91   2.10
## chr16    1.74    1.99   2.06   2.04     0.13    2.11   2.28
## chr17    1.09    1.22   1.27   1.26     0.08    1.31   1.42
## chr18    1.99    2.95   3.18   3.16     0.56    3.53   4.04
## chr19    0.84    0.95   1.01   1.00     0.09    1.05   1.21
## chr20    1.23    1.41   1.45   1.45     0.11    1.50   1.65
## chr21    0.83    0.96   0.99   0.99     0.07    1.03   1.12
## chr22    2.31    2.49   2.62   2.64     0.22    2.74   3.17
## chr23    0.52    0.62   0.64   0.64     0.04    0.66   0.72
## chr24    0.95    1.10   1.18   1.15     0.09    1.21   1.33
## chr25    1.45    1.53   1.56   1.58     0.08    1.63   1.75
## chr26    0.43    0.48   0.50   0.49     0.03    0.51   0.55
## chr27    9.34   12.10  14.47  13.96     2.48   15.66  18.96
## chr28    1.81    2.19   2.35   2.35     0.27    2.58   2.82
## chr29    0.39    0.46   0.48   0.48     0.05    0.51   0.57
## chr30    0.71    0.86   0.91   0.92     0.09    0.98   1.09
## chr31    1.18    1.55   1.81   1.82     0.33    2.06   2.40
## chrX     1.53    1.70   1.76   1.74     0.09    1.80   1.87
## chrAll 100.00  100.00 100.00 100.00     0.00  100.00 100.00
```

### Summary HTSeq:  
**Number of reads processed**


```r
processed.htseq <- rowSums(data.frame(sumHTSeq$Paired_Processed*2, 
    sumHTSeq$Single_Processed))
names(processed.htseq) <- animID
summSD(processed.htseq/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  32.2569  44.1879  49.2568  49.2636   7.2077  52.7932  63.8192
```

**Number of reads processed with no feature**


```r
no.feature.htseq <- rowSums(data.frame(sumHTSeq$Paired_No_Feature*2, 
    sumHTSeq$Single_No_Feature))
names(no.feature.htseq) <- animID
summSD(no.feature.htseq/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   5.9710   7.7088   8.6055   8.8993   1.6873  10.3501  13.2036
```

**Percent of reads processed with no feature atribute (no counts obtained) per animal**


```r
summSD(no.feature.htseq/processed.htseq*100)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  13.4120  16.0500  18.5109  18.1051   2.3407  19.7533  21.9180
```

**Percent of retained processed reads per animal**


```r
summSD(100 - (no.feature.htseq/processed.htseq*100))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  78.0820  80.2467  81.4891  81.8949   2.3407  83.9500  86.5880
```

### Summary expressed genes  
 Merge paired and single read counts per animal


```r
anim <- unlist(lapply(strsplit(colnames(paired), "_"), 
    function(x) x[[1]][1]))
counts <- do.call(cbind, lapply(1:ncol(paired), 
    function(x) rowSums(data.frame(paired[,x], single[,x]))))
colnames(counts) <- anim
rownames(counts) <- rownames(paired)
dim(counts)
```

```
## [1] 37875    23
```

**Remove HTSeq statistics from counts (laste five rows)**


```r
tail(counts[,1:5])
```

```
##                         X12401  X12402  X12403  X12610  X12611
## XLOC_037904               1026     887     678    1140    1207
## __no_feature           6811398 5479738 3852800 5187391 5972106
## __ambiguous                  0       0       0       0       0
## __too_low_aQual              0       0       0       0       0
## __not_aligned                0       0       0       0       0
## __alignment_not_unique       0       0       0       0       0
```

```r
counts <- counts[1:(nrow(counts)-5),]
tail(counts[,1:5])
```

```
##             X12401 X12402 X12403 X12610 X12611
## XLOC_037899      0      0      0      0      0
## XLOC_037900    229    278    193    240    324
## XLOC_037901      0      0      2      1      0
## XLOC_037902     28     51      6     20     62
## XLOC_037903      0      0      0      0      0
## XLOC_037904   1026    887    678   1140   1207
```

**Total number of expressed genes (genes with at least one count)**


```r
counts <- counts[rowSums(counts) > 0,]
nrow(counts)
```

```
## [1] 28625
```

**Number of expressed genes per chromosome**


```r
table(as.character(annot[rownames(counts), "chr"]))
```

```
## 
##  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 
##  2222  1485  1527   763   948   990   951   928   557   674   583  1538 
## chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29  chr3 chr30 
##   924   519   779   517   617   706   356   324   564   318  1224   318 
## chr31  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
##   240   986  1446  1252  1502  1123   735  1009
```

**Number of autosomal genes**


```r
sum(table(as.character(annot[rownames(counts), "chr"]))[-32])
```

```
## [1] 27616
```

### Summary genes for differential expression analysis   
**Number of gene transcripts with expression counts greater than 2 in all animals**  


```r
# Number of gene transcripts available for downstream analysis:
idx <- apply(counts,1, function(x) sum(x > 2) == ncol(counts))
counts <- counts[idx,]
nrow(counts)
```

```
## [1] 14155
```

Get annotation for expressed genes


```r
annot <- annot[rownames(counts),]
```

**Number of genes per chromosome**


```r
chr.count <- table(as.character(annot$chr))
chr.count
```

```
## 
##  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 
##  1077   729   809   393   548   505   469   505   252   321   291   787 
## chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29  chr3 chr30 
##   415   271   381   249   324   367   137   150   280   137   585   143 
## chr31  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
##   117   465   674   602   784   543   359   486
```

**Number of autosomal transcripts**


```r
sum(chr.count[-(length(chr.count))])
```

```
## [1] 13669
```

**Save count HTSeq table**


```r
write.table(counts, paste(dir, "HTSeq", "htseq_counts_RER_Thoroughbred.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Save RNA-Seq Pipeline read summary


```r
save(startPE, trimm.out, cond.out, tophat.out, uniq.out, depth, 
   depth.chr, ref, ratio.chr, processed.htseq, no.feature.htseq,
   counts, file=paste(getwd(), "retained_read_stats_RER_Thoroughbred.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Descriptive_Analysis_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RNA-Seq Pipeline Statistics
```

