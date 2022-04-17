## Meta-analysis
- **Project:** Assessment of Parkinson's disease pathogenic genes and risk variants
- **Author(s):** Cornelis, Mike, Andy, and Vanessa
- **Date Last Updated:** January 2022

### Quick Description: 
METAL is a tool for meta-analysis genomewide association scans. METAL can combine either (a) test statistics and standard errors or (b) p-values across studies (taking sample size and direction of effect into account). METAL analysis is a convenient alternative to a direct analysis of merged data from multiple studies. It is especially appropriate when data from the individual studies cannot be analyzed together because of differences in ethnicity, phenotype distribution, gender or constraints in sharing of individual level data imposed. Meta-analysis results in little or no loss of efficiency compared to analysis of a combined dataset including data from all individual studies.

## Workflow
0. Getting set up
1. Prepare files to match
2. Create METAL file
3. Run meta analysis
4. Annotate
5. File edits
6. Plot


## 0. Getting set up
```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
sinteractive --mem=100g --cpus-per-task=20
```

## 1. Prepare files to match

Add weights (number of people) to 23andme file
```
module load R
R
library(tidyverse)
library(data.table)
andme = fread("results_collaborators_all_963.txt")

# split alleles
andme = andme %>% separate("alleles", c("A1", "A2"), sep = "/")

# derive weight, i.e how many samples analyzed
andme =andme %>% dplyr::mutate(andme, "N_controls" = rowSums((andme[,12:14]), na.rm = TRUE))
andme =andme %>% dplyr::mutate(andme, "N_cases" = rowSums((andme[,17:19]), na.rm = TRUE))
andme %>% summarise_at(vars(N_controls, N_cases), funs(max), na.rm = T)
#  N_controls N_cases
#    3065473   25036

3065473+25034
[1] 3090507
andme$N_INFORMATIVE = 3090507

write.table(andme, "formeta_results_collaborators_all_963.txt", row.names= F, sep = "\t", quote =F)
```

### Edit AMP-PD x LNG to match 23andme layout
```
# start with AMP file
AMP = fread("AMP_PD_963_23andme_withcovars_score.SingleScore.assoc")

#edit layout
AMP$CHROM <- paste0("chr",AMP$CHROM)
AMP$REF1 =AMP$REF
AMP$ALT1 =AMP$ALT


AMP = AMP %>% unite("CHR.START.REF.ALT", c("CHROM", "POS", "REF", "ALT"), sep = ":")

#filter steps
#data = data %>% filter(Beta <5 & Beta > -5 & !is.na(Pvalue))

AMP$minorAllele <- ifelse(AMP$AF <= 0.5, as.character(AMP$ALT1), as.character(AMP$REF1))
AMP$majorAllele <- ifelse(AMP$AF <= 0.5, as.character(AMP$REF1), as.character(AMP$ALT1))

AMP = AMP %>% rename("ALT" = ALT1, "REF" = REF1)

write.table(AMP, "toMETA_SCORE_AMP.txt", quote = F, sep = "\t", row.names = F)
```

### Edit UKB to match 23andme layout
```
UKB = fread("UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc")

#edit layout
UKB$CHROM = paste0("chr",UKB$CHROM)
UKB$REF1 =UKB$REF
UKB$ALT1 =UKB$ALT

UKB = UKB %>% unite("CHR.START.REF.ALT", c("CHROM", "POS", "REF", "ALT"), sep = ":")

# filter steps
# data = data %>% filter(Beta <5 & Beta > -5 & !is.na(Pvalue))

UKB$minorAllele <- ifelse(UKB$AF <= 0.5, as.character(UKB$ALT1), as.character(AMP$REF1))
UKB$majorAllele <- ifelse(UKB$AF <= 0.5, as.character(UKB$REF1), as.character(AMP$ALT1))

UKB = UKB %>% rename("ALT" = ALT1, "REF" = REF1)

write.table(UKB, "toMETA_SCORE_UKBALL.txt", quote = F, sep = "\t", row.names = F)
```

### Edit 23andme
```
# make 23andme data match this layout
andme = fread("formeta_results_collaborators_all_963.txt")
andme$Alleles = andme$CHR.BP.REF.ALT
andme = andme %>% separate(Alleles, c("CHR", "BP", "REF", "ALT"), sep = ":")

andme = andme %>% select(-c(CHR, BP))

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

# Karl: I'd recommend keeping most of these thresholds, but experiment with relaxing p.date and p.batch. LRRK2 G2019S doesn't fail by much, so it may be sensible to use a slightly less stringent threshold there (instead of removing these filters altogether).

dim(andme)
[1] 869  40
```

Some imputed variants have the allele frequencies in a different column, so check this and make one frequency column
```
# overall how many are genotyped vs imputed?
andme %>% group_by(src) %>% tally()
# A tibble: 2 × 2
  src       n
  <chr> <int>
1 G       365
2 I       504

# correct those imputed ones with missing freq.b
andme$FREQ = andme$freq.b
andme_edit = andme %>% mutate(FREQ = if_else(src == "I" & is.na(freq.b), dose.b, FREQ))
andme_edit %>% select(CHR.BP.REF.ALT, src, dose.b, freq.b, FREQ)

write.table(andme, "toMETA_23andme_summary.txt", row.names = F, sep = "\t", quote=F)
```

Annotate all files for future use
```
R
library(dplyr)
library(tidyr)
ANDME = read.table("toMETA_23andme_summary.txt", header = T, sep = "\t")
AMP = read.table("toMETA_SCORE_AMP.txt", header = T, sep = "\t")
UKB = read.table("toMETA_SCORE_UKBALL.txt", header = T, sep = "\t")

ANDME_annotate = ANDME %>% select(CHR.BP.REF.ALT) %>% separate(CHR.BP.REF.ALT, c("CHR", "START", "REF", "ALT"), sep = ":")
ANDME_annotate$END = ANDME_annotate$START
ANDME_annotate= ANDME_annotate %>% select(CHR, START, END, REF, ALT)
names(ANDME_annotate) <- NULL
write.table(ANDME_annotate, "to_annotate_ANDME.txt", row.names = F, sep = "\t", quote =F) 

AMP_annotate = AMP %>% select(CHR.START.REF.ALT) %>% separate(CHR.START.REF.ALT, c("CHR", "START", "REF", "ALT"), sep = ":")
AMP_annotate$END = AMP_annotate$START
AMP_annotate= AMP_annotate %>% select(CHR, START, END, REF, ALT)
names(AMP_annotate) <- NULL
write.table(AMP_annotate, "to_annotate_AMP.txt", row.names = F, sep = "\t", quote =F) 

UKB_annotate = UKB %>% select(CHR.START.REF.ALT) %>% separate(CHR.START.REF.ALT, c("CHR", "START", "REF", "ALT"), sep = ":")
UKB_annotate$END = UKB_annotate$START
UKB_annotate= UKB_annotate %>% select(CHR, START, END, REF, ALT)
names(UKB_annotate) <- NULL
write.table(UKB_annotate, "to_annotate_UKB.txt", row.names = F, sep = "\t", quote =F) 
```

```
module load annovar

#gene build hg38

table_annovar.pl to_annotate_ANDME.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -outfile 23andme_annotated -nastring .

table_annovar.pl to_annotate_AMP.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -outfile AMP_annotated -nastring .

table_annovar.pl to_annotate_UKB.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -outfile UKB_annotated -nastring .

## new files
UKB_annotated.hg38_multianno.txt
AMP_annotated.hg38_multianno.txt
23andme_annotated.hg38_multianno.txt
```

Merge annotation back on files
```
R
library(tidyverse)

**23andme**
ANDME = read.table("toMETA_23andme_summary.txt", header = T, sep = "\t")
dim(ANDME)
# 869 41

ANDME_anno = read.table("23andme_annotated.hg38_multianno.txt", header = T, sep = "\t")
dim(ANDME_anno)
# 869 8
ANDME_anno = ANDME_anno %>% unite("CHR.BP.REF.ALT", c(Chr, Start, Ref, Alt), sep = ":")
ANDME_total = full_join(ANDME, ANDME_anno)
dim(ANDME_total)
# 869 48
write.table(ANDME_total, "23andme_summarystats_annotated.txt", row.names = F, sep = "\t", quote =F)

**AMP-PD**
AMP = read.table("toMETA_SCORE_AMP.txt", header = T, sep = "\t")
dim(AMP)
# 349 14

AMP_anno = read.table("AMP_annotated.hg38_multianno.txt", header = T, sep = "\t")
dim(AMP_anno)
# 349 11
AMP_anno = AMP_anno %>% unite("CHR.START.REF.ALT", c(Chr, Start, Ref, Alt), sep = ":")
AMP_total = full_join(AMP, AMP_anno)
dim(AMP_total)
# 349 21
write.table(AMP_total, "AMP_summarystats_annotated.txt", row.names = F, sep = "\t", quote =F)

**UKB**
UKB = read.table("toMETA_SCORE_UKBALL.txt", header = T, sep = "\t")
dim(UKB)
# 718 14

UKB_anno = read.table("UKB_annotated.hg38_multianno.txt", header = T, sep = "\t")
dim(UKB_anno)
# 718 11
UKB_anno = UKB_anno %>% unite("CHR.START.REF.ALT", c(Chr, Start, Ref, Alt), sep = ":")
UKB_total = full_join(UKB, UKB_anno)
dim(UKB_total)
# 718 21
write.table(UKB_total, "UKB_summarystats_annotated.txt", row.names = F, sep = "\t", quote =F)
```

Edit AAChange for individual plots
```
Follow script I wrote https://hackmd.io/pEyV63whT_m0KnGAWmPsnQ?view
# resulting files: Edited_AAChange_23andme_sumstats.txt,
Edited_AAChange_AMP_sumstats.txt,
Edited_AAChange_UKB_sumstats.txt
```

# 2. Create METAL file
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
MARKER CHR.START.REF.ALT
ALLELE minorAllele majorAllele
FREQ AF
EFFECT EFFECT
STDERR SE
PVALUE PVALUE
WEIGHT N_INFORMATIVE 
PROCESS toMETA_SCORE_AMP.txt

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER CHR.START.REF.ALT
ALLELE minorAllele majorAllele
FREQ AF
EFFECT EFFECT
STDERR SE
PVALUE PVALUE
WEIGHT N_INFORMATIVE 
PROCESS toMETA_SCORE_UKBALL.txt

# === DESCRIBE AND PROCESS THE THIRD INPUT FILE ===
MARKER CHR.BP.REF.ALT
ALLELE A2 A1
FREQ FREQ
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

# 3. Run meta-analysis
```
module load metal
metal my_METAL.txt

## Running second pass analysis to evaluate heterogeneity...
## Processing file 'toMETA_23andme_summary.txt'
## Processing file 'toMETA_SCORE_UKBALL.txt'
## Processing file 'toMETA_SCORE_AMP.txt'

## Executing meta-analysis ...
## Complete results will be stored in file 'MY_META_AMP_UKB_23andme1.tbl'
## Column descriptions will be stored in file 'MY_META_AMP_UKB_23andme1.tbl.info'
## Completed meta-analysis for 833 markers!
## Smallest p-value is 3.06e-384 at marker 'chr12:40340400:G:A'
```

Check output files and make plots
What do you need for that? (https://github.com/ipdgc/Manhattan-Plotter, https://github.com/GP2-TNC-WG/GP2-Bioinformatics-course/blob/master/Module_III.md#6-data-visualization-manhattan-plot-1)
- META analysis .tbl file
- ANNOTATE file (including Gene name and status)

# 4. Annotation 
### Prepare file for annotation
```
module load R
R
library(tidyverse)
meta = read.table("MY_META_AMP_UKB_23andme1.tbl", header = T, sep = "\t")

#create file for ANNOTATION
annotate = meta %>% select(MarkerName)
annotate = annotate %>% separate("MarkerName", c("CHR", "START", "REF", "ALT"), sep = ":")
annotate$END = annotate$START

annotate = annotate %>% select(CHR, START, END, REF, ALT)
#remove header
names(annotate)<-NULL

write.table(annotate, "META_to_annotate.txt", quote = F, sep = "\t", row.names = F)
q()
n
```

### Run ANNOVAR
```
module load annovar

#gene build hg38

table_annovar.pl META_to_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -outfile META_rare_variants_first_run -nastring .

## new files
META_rare_variants_first_run.hg38_multianno.txt
```

# 5. File edits
### Add OR and 95% CI
```
R
library(tidyverse)

Genename = read.table("META_rare_variants_first_run.hg38_multianno.txt", header = T, sep = "\t")
meta = read.table("MY_META_AMP_UKB_23andme1.tbl", header = T, sep = "\t")


Genename = Genename %>% unite("CHR.BP.REF.ALT", c("Chr", "Start", "Ref", "Alt"), sep = ":")
Genename = Genename %>% select(CHR.BP.REF.ALT, Gene.refGene, AAChange.refGene)
Genename = Genename %>% rename("Gene" = Gene.refGene, "AAChange" = AAChange.refGene)

meta = meta %>% rename("CHR.BP.REF.ALT" = MarkerName)

leftjoin = left_join(meta, Genename)
dim(leftjoin)
# 833 17

write.table(leftjoin, "META_AMP_UKB_23andme.txt", row.names =F, sep ="\t", quote = F)
```

### Edit AAChange 
```
Follow script I wrote https://hackmd.io/pEyV63whT_m0KnGAWmPsnQ?view
# resulting file: Edited_AAChange_Metaanalysis.txt
```

# 6. Plot data and write result files
```
library(ggplot2)
library(data.table)
library(tidyverse)
library(forcats)
library(stringr)

meta = fread("Edited_AAChange_Metaanalysis.txt")

# flip negative Effects (betas)
meta = meta %>% mutate(Effectv2 = if_else(Effect <0, abs(Effect), Effect))

# add OR and 95% CI
meta = meta %>% mutate(OR = exp(Effectv2), L95 = exp(Effectv2 - 1.96*StdErr), U95 = exp(Effectv2 + 1.96*StdErr))
meta$U95 = as.numeric(meta$U95)
```

Write file of all variants passing p<0.05
```
meta0.05 = meta%>% filter(P.value <0.05)
meta0.05 = tidyr::separate(meta0.05, "MarkerName", c("CHR", "BP", "REF", "ALT"), sep = ":")
write.table(meta0.05, "Results_meta_analysis_pvalue0.05.txt", row.names = F, sep = "\t", quote = F)

# add width of CI to filter those that are way too large for plotting
meta0.05 = meta0.05 %>% mutate(CI_range = U95-L95)
CI_filter = meta0.05 %>% filter(CI_range<30)

# edit those with multiple AAChanges - choose smallest NM_transcript
CI_filter = CI_filter %>% mutate(VariantName = if_else(VariantName == "FBXO7_M36I, FBXO7_M115I", "FBXO7_M115I",
                                                           if_else(VariantName == "DNAJC13_A1463S, DNAJC13_A1468S", "DNAJC13_A1463S",
                                                                   if_else(VariantName == "EIF4G1_M432V, EIF4G1_M439V", "EIF4G1_M432V",
                                                                           if_else(VariantName == "EIF4G1_I807V, EIF4G1_I813V, EIF4G1_I806V", "EIF4G1_I806V",
                                                                                   if_else(VariantName == "EIF4G1_L604V, EIF4G1_L611V", "EIF4G1_L604V",
                                                                                           if_else(VariantName == "DNAJC6_V353M, DNAJC6_V296M, DNAJC6_V283M", "DNAJC13_V296M",
                                                                                                   if_else(VariantName == "DNAJC13_R1830H, DNAJC13_R1835H", "RNAJC13_R1830H", VariantName))))))))


# plot
CI_filtered = CI_filter %>% 
  ggplot(CI_filter, mapping = aes(x= OR, y = reorder(VariantName, -OR)))+
  geom_vline(aes(xintercept =1), size = .5, linetype = "dashed")+
  geom_errorbarh(aes(xmax = U95, xmin = L95), size = .5, height = .2) +
  geom_point(size = 3.5, aes(color = OR)) +
  scale_x_continuous(breaks = seq(0,25,1.5), labels = seq(0,25,1.5), limits = c(0,25)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Variants")+
  xlab("OR (95% CI)")+
  ggtitle("Meta-analysis estimates p<0.05")+
  theme(plot.title = element_text(hjust=0.5))
  
ggsave("Results_ForestPlot_23andme_META_p0.05.png", CI_filtered, width = 12, height = 10, dpi=300, units = "in")
```

Write file with those with smalles p-values
```
## filtering p-values, otherwise too many to plot
meta$log10Praw <- -1*log(meta$P.value, base = 10)
meta$log10P <- ifelse(meta$log10Praw>40, 40, meta$log10Praw) 
gwasFiltered <- subset(meta, log10P > 3.114074)
gwasFiltered %>% summarise(min = min(OR),
                           max = max(OR))

selected = gwasFiltered 
write.table(selected, "Results_meta_analysis_selectedpvalue.txt", row.names = F, sep = "\t", quote = F)

# remove wide CIs - visually
gwasFiltered = gwasFiltered %>% filter(VariantName != "LRRK2_I2020T" & VariantName != "GBA_S146L" & VariantName != "FBXO7_T22M" & VariantName != "PRKN_R33X" & VariantName != "LRRK2_R1441H")

filtered = gwasFiltered %>% 
  ggplot(gwasFiltered, mapping = aes(x= OR, y = reorder(VariantName, -OR)))+
  geom_vline(aes(xintercept =1), size = .5, linetype = "dashed")+
  geom_errorbarh(aes(xmax = U95, xmin = L95), size = .5, height = .2) +
  geom_point(size = 3.5, aes(color = OR)) +
  scale_x_continuous(breaks = seq(0,14,1.5), labels = seq(0,14,1.5), limits = c(0,14)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Variants")+
  xlab("OR (95% CI)")+
  ggtitle("Meta-analysis estimates of most significant variants")+
  theme(plot.title = element_text(hjust=0.5))

ggsave("Results_ForestPlot_23andme_META.png", filtered, width = 12, height = 5, dpi=300, units = "in")
```
