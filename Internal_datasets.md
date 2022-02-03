## Extract 1091 variants from internal data sets
- Vanessa

---
### Quick Description: 
Extracting 963 variants (available in GRCH38 from 23andme data) from internal data sets and generate association stats on those variants only.
Follow Mary's script for burden testing but instead of burden testing, do --score test

Extracted variants from 23andme in GRCh37, AMP-PD are in GRCh38

Connect to biowulf and start node
```
sinteractive --mem=100g --cpus-per-task=20
# new folder
/data/CARD/projects/23andme_annotation/variants_in_internal_datasets
```

Create new variant file from old and new 23and me
```
module load R
R
library(dplyr)
library(tidyr)

failed = read.table("results_for_collaborators_failed_qc.csv", header = T, sep = ",")
dim(failed)
# 459 33

old = read.table("results_for_collaborators.csv", sep = ",", header = T)
dim(old)
# 504 33

fulljoin_old_failed = full_join(failed, old) # looks like they're all new
dim(fulljoin_old_failed)
# 963 33
head(fulljoin_old_failed)
fulljoin_old_failed %>% group_by(position) %>% tally() %>% filter(n>1) #no duplicates

write.table(fulljoin_old_failed, "results_collaborators_all_963.txt", row.names = F, quote = F, sep = "\t")

# create file with only IDs, so I can grep them from data sets
fulljoin_old_failed = fulljoin_old_failed %>% select(scaffold, position)
fulljoin_old_failed = fulljoin_old_failed %>% unite(combo, c("scaffold", "position"), sep = ":")
names(fulljoin_old_failed) = NULL
head(fulljoin_old_failed)

write.table(fulljoin_old_failed, "963_variants_in23andme_variants_IDs_only.txt", row.names = F, quote = F, sep = "\t")
q()
n
```
All files now in GRCH38!

## Call on data sets and extract variants

### AMP-PD 
```
cd /data/CARD/PD/AMP_NIH/no_relateds/
# /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins*

#covariate file for burdens
/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt
```

Grep variants and write into new text file
```
grep -w -f 963_variants_in23andme_variants_IDs_only.txt /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.bim > List_963variants_alleles_AMPPD.txt

wc -l List_963variants_alleles_AMPPD.txt
# 367 List_963variants_alleles_AMPPD.txt
```

Double check if these are the correct variants
```
R
library(dplyr)
library(tidyr)
variants_new = read.table("List_963variants_alleles_AMPPD.txt", header =F, sep = "\t")
variants = read.table("963_variants_in23andme_variants_IDs_only.txt",  header =F, sep = "\t")

variants_new = variants_new %>% separate(V2, c("chr", "position", "allele1", "allele2"), sep = ":")
Check = variants_new %>% select(chr, position)
Check = Check %>% unite("V1",c(chr:position), sep = ":")
join = left_join(Check, variants)
dim(join)
# 367 1
q()
n
```

Write new binary files - only 19 variants
```
module load plink

plink --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins --extract List_963variants_alleles_AMPPD.txt --make-bed --out AMP_PD_963_23andme

###
257652 MB RAM detected; reserving 128826 MB for main workspace.
47899791 variants loaded from .bim file.
7986 people (4345 males, 3641 females) loaded from .fam.
7986 phenotype values loaded from .fam.
--extract: 367 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 7986 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999999.
367 variants and 7986 people pass filters and QC.
Among remaining phenotypes, 3376 are cases and 4610 are controls.
--make-bed to AMP_PD_963_23andme.bed + AMP_PD_963_23andme.bim +
AMP_PD_963_23andme.fam ... done.
```
AMP-PD only contains 367 out of 963 variants provided by 23andme

Generate frequency file
```
module load plink
plink --bfile AMP_PD_963_23andme --freq --out AMP_PD_963_23andme

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
367 variants loaded from .bim file.
7986 people (4345 males, 3641 females) loaded from .fam.
7986 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 7986 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999999.
--freq: Allele frequencies (founders only) written to AMP_PD_963_23andme.frq .
```
Run association test in rvtest
```
# convert plink files to vcf for rvtest
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile AMP_PD_963_23andme \
--export vcf bgz id-paste=iid --out AMP_PD_963_23andme${OUTNAME} --mac 1

tabix -p vcf  AMP_PD_963_23andme${OUTNAME}.vcf.gz

# new files
AMP_PD_963_23andme.vcf.gz
AMP_PD_963_23andme.vcf.gz.tbi
```

Run rvtests
```
module load rvtests
rvtest --inVcf AMP_PD_963_23andme.vcf.gz --pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --pheno-name PD_PHENO --out AMP_PD_963_23andme_withcovars_score --single wald,score --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5

# new files
AMP_PD_963_23andme_withcovars_score.log
AMP_PD_963_23andme_withcovars_score.SingleScore.assoc
AMP_PD_963_23andme_withcovars_score.SingleWald.assoc
```

### UKB
```
# PLINK file merged with no relateds
/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass.*
# these files are in plink2 format - convert

module load plink/2
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass --make-bed --out MERGED_UKB_first_pass

grep -w -f 963_variants_in23andme_variants_IDs_only.txt MERGED_UKB_first_pass.bim > List_963variants_alleles_UKB.txt
wc -l List_963variants_alleles_UKB.txt
# 843 List_963variants_alleles_UKB.txt
```

Double check if these are the correct variants
```
R
library(dplyr)
library(tidyr)
variants_new = read.table("List_963variants_alleles_UKB.txt", header =F, sep = "\t")
variants = read.table("List_1091variants_ID_only.txt",  header =F, sep = "\t")

variants_new = variants_new %>% separate(V2, c("chr", "position", "allele1", "allele2"), sep = ":")
Check = variants_new %>% dplyr::select(chr, position)
Check = Check %>% unite("V1",c(chr:position), sep = ":")
join = left_join(Check, variants)
dim(join)
# 843 1
q()
n
```

Write new binary files - 843 variants

```
module load plink/1.9.0-beta4.4

plink --bfile MERGED_UKB_first_pass --extract List_963variants_alleles_UKB.txt --make-bed --out UKB_963_23andme

###
257652 MB RAM detected; reserving 128826 MB for main workspace.
16285684 variants loaded from .bim file.
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to UKB_963_23andme.nosex .
--extract: 843 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
843 variants and 200648 people pass filters and QC.
Note: No phenotypes present.
--make-bed to UKB_963_23andme.bed + UKB_963_23andme.bim + UKB_963_23andme.fam
... done.
```

Generate frequency file
```
module load plink
plink --bfile UKB_963_23andme --freq --out UKB_963_23andme

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
843 variants loaded from .bim file.
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to UKB_963_23andme.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
--freq: Allele frequencies (founders only) written to UKB_963_23andme.frq .

```

Run association test in rvtest
```
# convert binary files to vcf for input in rvtest
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile UKB_963_23andme \
--export vcf bgz id-paste=iid --out UKB_963_23andme${OUTNAME} --mac 1

tabix -p vcf  UKB_963_23andme${OUTNAME}.vcf.gz

# new files written: 
UKB_963_23andme.vcf.gz
UKB_963_23andme.vcf.gz.tbi
```

Run rvtests using different phenotype and covariate files (4 different)
```
module load rvtests

# using UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.t
rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleWald.assoc


# using UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleWald.assoc


# using UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files 
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleWald.assoc


# using UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt
rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_SIBLING_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_PD_SIBLING_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_SIBLING_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_SIBLING_CONTROL_2021_with_PC.SingleWald.assoc
```

## Merge files together

```
R
library(dplyr)
library(tidyr)


#####
andme = read.table("results_collaborators_all_963.txt", sep = "\t", header = T)

# layout
| assay.name  | scaffold | position  | alleles |
|-------------|----------|-----------|---------|
| rs536425950 | chr1     | 155235024 | C/T     |

# change
andme = andme %>% rename("Chr" = scaffold, "Start" = position)
andme$End = andme$Start
andme = andme %>% separate(alleles, c("REF", "ALT"), sep = "/")
andme = andme %>% select(Chr, Start, End, REF, ALT, assay.name, pvalue:p.batch)
dim(andme)
# 963  35

#####
AMP_score = read.table("AMP_PD_963_23andme_withcovars_score.SingleScore.assoc", sep = "\t", header = T)
# layout
| CHROM | POS     | REF | ALT | N_INFORMATIVE |
|-------|---------|-----|-----|---------------|
| 1     | 7965399 | G   | A   | 7974          |

# change
AMP_score = AMP_score %>% rename("Chr" = CHROM, "Start" = POS)
AMP_score$End = AMP_score$Start
AMP_score$Chr = paste0('chr', AMP_score$Chr)
AMP_score = AMP_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
AMP_score = AMP_score %>% rename_with(~paste0(., "_AMP_score"), N_INFORMATIVE:PVALUE)
dim(AMP_score)
# 351  14

# layout we need for annotation
| Chr  | Start     | End       | Ref | Alt |
|------|-----------|-----------|-----|-----|
| chr1 | 155235006 | 155235006 | A   | G   |


#####
UKB_ALL_score = read.table("UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc",  sep = "\t", header = T)

# change
UKB_ALL_score = UKB_ALL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_ALL_score$End = UKB_ALL_score$Start
UKB_ALL_score$Chr = paste0('chr', UKB_ALL_score$Chr)
UKB_ALL_score = UKB_ALL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_ALL_score = UKB_ALL_score %>% rename_with(~paste0(., "_UKB_ALL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_ALL_score)
# 843  14

#####
UKB_CASE_CTRL_score = read.table("UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleScore.assoc", sep = "\t", header = T)

# change
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_CASE_CTRL_score$End = UKB_CASE_CTRL_score$Start
UKB_CASE_CTRL_score$Chr = paste0('chr', UKB_CASE_CTRL_score$Chr)
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% rename_with(~paste0(., "_UKB_CASE_CONTROL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_CASE_CTRL_score)
# 843  14

#####
UKB_PARENT_CTRL_score = read.table("UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleScore.assoc", sep = "\t", header = T)

# change
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_PARENT_CTRL_score$End = UKB_PARENT_CTRL_score$Start
UKB_PARENT_CTRL_score$Chr = paste0('chr', UKB_PARENT_CTRL_score$Chr)
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% rename_with(~paste0(., "_UKB_PARENT_CONTROL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_PARENT_CTRL_score)
# 843  14


#####
UKB_SIBLING_CTR_score = read.table("UKB_963_23andme_withcovars_score_PD_SIBLING_CONTROL_2021_with_PC.SingleScore.assoc", sep = "\t", header = T)

# change
UKB_SIBLING_CTR_score = UKB_SIBLING_CTR_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_SIBLING_CTR_score$End = UKB_SIBLING_CTR_score$Start
UKB_SIBLING_CTR_score$Chr = paste0('chr', UKB_SIBLING_CTR_score$Chr)
UKB_SIBLING_CTR_score = UKB_SIBLING_CTR_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_SIBLING_CTR_score = UKB_SIBLING_CTR_score %>% rename_with(~paste0(., "_UKB_SIBLING_CONTROL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_SIBLING_CTR_score)
# 843  14

```


Start merging
```
merge1 = left_join(andme, AMP_score)
dim(merge1)
# 963 44

merge2 = left_join(merge1, UKB_ALL_score)
dim(merge2)
# 963 53

merge3 = left_join(merge2, UKB_CASE_CTRL_score)
dim(merge3)
# 963 62

merge4 = left_join(merge3, UKB_PARENT_CTRL_score)
dim(merge4)
# 963 71

merge5 = left_join(merge4, UKB_SIBLING_CTR_score)
dim(merge5)
# 963 80

# check for merging errors, i.e. are any columns all NA
sapply(merge5, function(x)all(is.na(x)))
# no, all good

write.table(merge5, "963variants_23andme_AMPPD_UKB_cov_merged_score.txt", row.names = F, sep = "\t", quote= F)
```

Write variant file for annotation
```
to_annotate = merge5 %>% select(Chr, Start, End, REF, ALT)
names(to_annotate) = NULL
write.table(to_annotate, "963variants_to_annotate.txt", sep = "\t", row.names =F, quote= F)
```

## Annotate 
```
module load annovar

#gene build hg38

table_annovar.pl 963variants_to_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150 \
-operation g,f -outfile 963variants_23_AMP_UKB -nastring .

## new files
963variants_23_AMP_UKB.hg38_multianno.txt
```

## Final file
Add annotation (Gene name) to merged file

```
R
library(dplyr)
library(tidyr)

variants = read.table("963variants_23andme_AMPPD_UKB_cov_merged_score.txt", header = T, sep = "\t")

Genenames = read.table("963variants_23_AMP_UKB.hg38_multianno.txt", header = T, sep = "\t")
Genenames = Genenames %>% select(Start, Gene.refGene)

#merge Gene name onto list
mergename = left_join(variants, Genenames)
```

Add frequencies
```
AMP_freq = read.table("AMP_PD_963_23andme.frq", header = T, sep = "")
AMP_freq = AMP_freq %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
AMP_freq = AMP_freq %>% select(Start, MAF) %>% rename("AMP_MAF" = MAF) 

mergename$Start = as.character(mergename$Start)
mergename_freq1 = left_join(mergename, AMP_freq)

UKB_freq = read.table("UKB_963_23andme.frq", header = T, sep = "")
UKB_freq = UKB_freq %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
UKB_freq = UKB_freq %>% select(Start, MAF) %>% rename("UKB_MAF" = MAF) 
mergename_freq1_2 = left_join(mergename_freq1, UKB_freq)

# re-arrange file a bi
mergename_freq1_2 = mergename_freq1_2 %>% select(Chr:assay.name, Gene.refGene, pvalue:p.batch, AMP_MAF, N_INFORMATIVE_AMP_score:PVALUE_AMP_score, UKB_MAF, N_INFORMATIVE_UKB_ALL_score:PVALUE_UKB_SIBLING_CONTROL_score) 

write.table(mergename_freq1_2, "963_variants_score_AMP_UKB_wide.txt", row.names = F, sep = "\t", quote = F)

```


