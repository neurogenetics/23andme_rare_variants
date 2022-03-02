## Meta-analysis
- **Project:** Assessment of Parkinson's disease pathogenic genes and risk variants
- **Author(s):** Cornelis, Mike, Andy, and Vanessa
- **Date Last Updated:** January 2022

### Quick Description: 
METAL is a tool for meta-analysis genomewide association scans. METAL can combine either (a) test statistics and standard errors or (b) p-values across studies (taking sample size and direction of effect into account). METAL analysis is a convenient alternative to a direct analysis of merged data from multiple studies. It is especially appropriate when data from the individual studies cannot be analyzed together because of differences in ethnicity, phenotype distribution, gender or constraints in sharing of individual level data imposed. Meta-analysis results in little or no loss of efficiency compared to analysis of a combined dataset including data from all individual studies.

### Working directory
```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
sinteractive --mem=100g --cpus-per-task=20
```

### Prepare files

Add weights (number of people) to 23andme file
```
module load R
R
library(tidyverse)
library(data.table)
andme = fread("results_collaborators_all_963.txt")

# derive weight, i.e how many samples analyzed
andme =andme %>% dplyr::mutate(andme, "N_controls" = rowSums((andme[,13:15]), na.rm = TRUE))
andme =andme %>% dplyr::mutate(andme, "N_cases" = rowSums((andme[,18:20]), na.rm = TRUE))
andme %>% summarise_at(vars(N_controls, N_cases), funs(max), na.rm = T)
#  N_controls N_cases
#    3065473   25034

3065473+25034
[1] 3090507
andme$N_INFORMATIVE = 3090507

andme$markerID <- paste(andme$scaffold,andme$position, sep = ":")

write.table(andme, "results_collaborators_all_963.txt", row.names= F, sep = "\t", quote =F)
```

Edit layout SCORE file and merge with info from 23andme
```
# start with AMP file
AMP = fread("AMP_PD_963_23andme_withcovars_score.SingleScore.assoc")

#edit layout
AMP$chr <- paste("chr",AMP$CHROM, sep = "")
AMP$markerID <- paste(AMP$chr,AMP$POS, sep = ":")
AMP = AMP %>% select(-c(CHROM, chr)) 

#filter steps
#data = data %>% filter(Beta <5 & Beta > -5 & !is.na(Pvalue))

AMP$minorAllele <- ifelse(AMP$AF <= 0.5, as.character(AMP$ALT), as.character(AMP$REF))
AMP$majorAllele <- ifelse(AMP$AF <= 0.5, as.character(AMP$REF), as.character(AMP$ALT))

#data$beta <- ifelse(data$freq.b <= 0.5, data$Beta, data$Beta*-1)
#data$se <- data$SE
#data$maf <- ifelse(data$freq.b <= 0.5, data$freq.b, 1 - data$freq.b)
#data$P <- data$Pvalue
#dat0 <- data[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]

write.table(AMP, "toMETA_SCORE_AMP.txt", quote = F, sep = "\t", row.names = F)
```

Repeat for UKB
```
UKB = fread("UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc")

#edit layout
UKB$chr <- paste("chr",UKB$CHROM, sep = "")
UKB$markerID <- paste(UKB$chr,UKB$POS, sep = ":")
UKB = UKB %>% select(-c(CHROM, chr)) 

#filter steps
#data = data %>% filter(Beta <5 & Beta > -5 & !is.na(Pvalue))

UKB$minorAllele <- ifelse(UKB$AF <= 0.5, as.character(UKB$ALT), as.character(AMP$REF))
UKB$majorAllele <- ifelse(UKB$AF <= 0.5, as.character(UKB$REF), as.character(AMP$ALT))

#data$beta <- ifelse(data$freq.b <= 0.5, data$Beta, data$Beta*-1)
#data$se <- data$SE
#data$maf <- ifelse(data$freq.b <= 0.5, data$freq.b, 1 - data$freq.b)
#data$P <- data$Pvalue
#dat0 <- data[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]

write.table(UKB, "toMETA_SCORE_UKBALL.txt", quote = F, sep = "\t", row.names = F)
```

Edit 23andme
```
# make 23andme data match this layout
andme = andme %>% select(markerID, assay.name, scaffold:N_INFORMATIVE)

# apply suggested thresholds (suggested by Karl Heilbron @23andme)
# this shows how many variants would be removed

## gt.rate < 0.9 (genotyping rate)
andme %>% filter(src == "G" & gt.rate <0.9) %>% tally()
#  n
#  5

andme %>% filter(src == "G" & p.date <1e-50) %>% tally()
#   n
#  14

andme %>% filter(src == "I" & avg.rsqr < 0.5) %>% tally()
#  n
#  91

andme %>% filter(src == "I" & min.rsqr < 0.5) %>% tally()
#   n
# 145

andme %>% filter(src == "I" & p.batch < 1e-50) %>% tally()
#  n
# 13

dim(andme)
[1] 963  38

# now apply them 
andme = andme %>% filter(src == "G" & gt.rate >0.9 | src == "G" & p.date >1e-50 | src == "I" & avg.rsqr > 0.5| src == "I" & min.rsqr > 0.5|src == "I" & p.batch > 1e-50)

dim(andme)
[1] 869  38

write.table(andme, "toMETA_23andme_summary.txt", row.names = F, sep = "\t", quote=F)

# Karl: I'd recommend keeping most of these thresholds, but experiment with relaxing p.date and p.batch. LRRK2 G2019S doesn't fail by much, so it may be sensible to use a slightly less stringent threshold there (instead of removing these filters altogether).
```

### Create metal file
Adapted from: https://github.com/neurogenetics/GWAS-pipeline
Create file and call it my_METAL.txt
```
#../generic-metal/metal metalAll.txt
#THIS SCRIPT EXECUTES AN ANALYSIS OF THREE STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile3.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
LABEL TotalSampleSize as N # If input files have a column for the sample size labeled as 'N'
# LOAD THE FIRST TWO INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ AF
EFFECT EFFECT
STDERR SE
PVALUE PVALUE
WEIGHT N_INFORMATIVE 
PROCESS toMETA_SCORE_AMP.txt

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ AF
EFFECT EFFECT
STDERR SE
PVALUE PVALUE
WEIGHT N_INFORMATIVE 
PROCESS toMETA_SCORE_UKBALL.txt

# === DESCRIBE AND PROCESS THE THIRD INPUT FILE ===
MARKER markerID
ALLELE A1 A2
FREQ freq.b
EFFECT effect
STDERR stderr
PVALUE pvalue
WEIGHT N_INFORMATIVE 
PROCESS toMETA_23andme_summary.txt

OUTFILE MY_META_AMP_UKB_23andme .tbl
MINWEIGHT 10000
ANALYZE HETEROGENEITY

QUIT
```

Then load module and run
```
module load metal
metal my_METAL.txt

## Executing meta-analysis ...
## Complete results will be stored in file 'MY_META_AMP_UKB_23andme1.tbl'
## Column descriptions will be stored in file 'MY_META_AMP_UKB_23andme1.tbl.info'
## Completed meta-analysis for 834 markers!
## Smallest p-value is 1.52e-315 at marker 'chr12:40340400'


Check output files and make plots
```
module load R
R
library(tidyverse)
library(data.table)

meta = fread("MY_META_AMP_UKB_23andme1.tbl")

meta$position = meta$MarkerName
meta = meta %>% separate(position, c("chr", "bp"), sep = ":")

PlotFile = meta %>% select(MarkerName, chr

SNP	CHR	BP	P

```

