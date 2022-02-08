## Extract 963 variants from internal data sets
- Vanessa

---
### Quick Description: 
Extracting 963 variants (available in GRCH38 from 23andme data) from internal data sets and generate association stats on those variants only.
Follow Mary's script for burden testing but instead of burden testing, do --score test
Files can also be found here: https://drive.google.com/drive/u/1/folders/18kISsSBgz7iD9Gu9000tQupiY-SD1TOH
Extracted variants from 23andme in GRCh37, AMP-PD are in GRCh38

Connect to biowulf and start node
```
sinteractive --mem=100g --cpus-per-task=20
# new folder
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
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

old = read.table("/data/CARD/projects/23andme_annotation/results_for_collaborators.csv", sep = ",", header = T)
dim(old)
# 504 33

fulljoin_old_failed = full_join(failed, old) 
# looks like they're all new (distinct)
dim(fulljoin_old_failed)
# 963 33
head(fulljoin_old_failed)
fulljoin_old_failed %>% group_by(position) %>% tally() %>% filter(n>1) #no duplicates

variants_1091 = read.table("List_1091variants_ID_only.txt", header = F, sep = "\t")

write.table(fulljoin_old_failed, "results_collaborators_all_963.txt", row.names = F, quote = F, sep = "\t")

# create file with only IDs, so I can grep them from data sets
fulljoin_old_failed = fulljoin_old_failed %>% select(scaffold, position)
fulljoin_old_failed = fulljoin_old_failed %>% unite(combo, c("scaffold", "position"), sep = ":")
names(fulljoin_old_failed) = "V1"
head(fulljoin_old_failed)

# check if these 23andme provided variants are all on our list
antijoin = anti_join(fulljoin_old_failed, variants_1091)
dim(antijoin)
[1] 963   1

names(antijoin) = NULL
write.table(antijoin, "963_variants_in23andme_variants_IDs_only.txt", row.names = F, quote = F, sep = "\t")
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
Frequency file for cases and controls
```
plink --bfile AMP_PD_963_23andme --assoc --pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --pheno-name PD_PHENO --out AMP_PD_963_23andme_pheno

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
367 variants loaded from .bim file.
7986 people (4345 males, 3641 females) loaded from .fam.
7986 phenotype values present after --pheno.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 7986 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999999.
367 variants and 7986 people pass filters and QC.
Among remaining phenotypes, 3376 are cases and 4610 are controls.
Writing C/C --assoc report to AMP_PD_963_23andme_pheno.assoc ...
done.
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

## .psam IDs
#IID	SEX
-000001	NA
-000002	NA

##
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

This file needs to be fixed because of an odd format where column 1 is all zeros (column 1 & 2 need to be same)
```
cut -d " " -f 2 UKB_963_23andme.fam > column2.txt
cut -d " " -f 2,3,4,5,6 UKB_963_23andme.fam > column23456.txt
scp UKB_963_23andme.fam UKB_963_23andme_ORIGINAL.fam
paste column2.txt column23456.txt > UKB_963_23andme.fam
```


Generate frequency file
```
module load plink
plink --bfile UKB_963_23andme --freq --out UKB_963_23andme

##
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to UKB_963_23andme.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
--freq: Allele frequencies (founders only) written to UKB_963_23andme.frq .
```

Frequency file for different phenotype groups
```
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
843 variants loaded from .bim file.
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.nosex .
45857 phenotype values present after --pheno.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
843 variants and 200648 people pass filters and QC.
Among remaining phenotypes, 7806 are cases and 38051 are controls.  (154791
phenotypes are missing.)
Writing C/C --assoc report to
UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.assoc ...
done.


##### PD_CASE_CONTROL
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
843 variants loaded from .bim file.
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC.nosex
.
6748 phenotype values present after --pheno.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
843 variants and 200648 people pass filters and QC.
Among remaining phenotypes, 1105 are cases and 5643 are controls.  (193900
phenotypes are missing.)
Writing C/C --assoc report to
UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC.assoc ...
done.


##### PD_PARENT_CONTROL_2021_with_PC
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC

##
386449 MB RAM detected; reserving 193224 MB for main workspace.
843 variants loaded from .bim file.
200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC.nosex .
34978 phenotype values present after --pheno.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 200648 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999133.
843 variants and 200648 people pass filters and QC.
Among remaining phenotypes, 6033 are cases and 28945 are controls.  (165670
phenotypes are missing.)
Writing C/C --assoc report to
UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC.assoc ...
done.
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

# using UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
# # A tibble: 4 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL    38051
2 parent      6033
3 PD          1105
4 sibling      668

rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleWald.assoc


# using UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
# # A tibble: 2 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL     5643
2 PD          1105

rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleWald.assoc


# using UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
# # A tibble: 2 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL    28945
2 parent      6033

rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files 
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleWald.assoc
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

# check for merging errors, i.e. are any columns all NA
sapply(merge4, function(x)all(is.na(x)))
# all good

write.table(merge4, "963variants_23andme_AMPPD_UKB_cov_merged_score.txt", row.names = F, sep = "\t", quote= F)
```

Write variant file for annotation
```
to_annotate = merge4 %>% select(Chr, Start, End, REF, ALT)
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
Genenames = Genenames %>% select(Start, Gene.refGene, AAChange.refGene)

#merge Gene name onto list
mergename = left_join(variants, Genenames)
dim(mergename)
# 963 73
```

Add AMP frequencies
```
# MAF for AMP
AMP_freq = read.table("AMP_PD_963_23andme.frq", header = T, sep = "")
| Chr | SNP              | A1 | A2 | MAF       | NCHROBS |
|-----|------------------|----|----|-----------|---------|
| 1   | chr1:7965399:G:A | A  | G  | 0.0001878 | 15972   |

## change layout
AMP_freq = AMP_freq %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
AMP_freq = AMP_freq %>% select(Start, MAF, A1, A2) %>% rename("AMP_MAF" = MAF)

mergename$Start = as.character(mergename$Start)
mergename_freq1 = left_join(mergename, AMP_freq)
mergename_freq1 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq1)
# 964 74

# MAF for phenotypes
AMP_phenoMAF = read.table("AMP_PD_963_23andme_pheno.assoc", header = T, sep = "")
| CHR | SNP              | BP      | A1 | F_A       | F_U       | A2 | CHISQ     | P      | OR     |
|-----|------------------|---------|----|-----------|-----------|----|-----------|--------|--------|
| 1   | chr1:7965399:G:A | 7965399 | A  | 0.0001481 | 0.0002169 | G  | 9.829e-02 | 0.7539 | 0.6827 |

AMP_phenoMAF = AMP_phenoMAF %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
AMP_phenoMAF$End = AMP_phenoMAF$Start
AMP_phenoMAF = AMP_phenoMAF %>% select(chr, Start, End, A1, A2, F_A, F_U, CHISQ, P, OR)
AMP_phenoMAF = AMP_phenoMAF %>% rename_with(~paste0(., "_AMP_PhenoFreq"), F_A:OR)
| chr  | Start   | End     | A1 | A2 | F_A_AMP_PhenoFreq | F_U_AMP_PhenoFreq | CHISQ_AMP_PhenoFreq | P_AMP_PhenoFreq | OR_AMP_PhenoFreq |
|------|---------|---------|----|----|-------------------|-------------------|---------------------|-----------------|------------------|
| chr1 | 7965399 | 7965399 | A  | G  | 0.0001481         | 0.0002169         | 9.829e-02           | 0.7539          | 0.6827           |

mergename_freq1$End = as.character(mergename_freq1$End)
mergename_freq2 = left_join(mergename_freq1, AMP_phenoMAF)

mergename_freq2 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq2)
# 966 81
```

Add UKB frequency (phenotype and MAF) files
```
UKB_freq = read.table("UKB_963_23andme.frq", header = T, sep = "")
| CHR | SNP              | A1 | A2 | MAF       | NCHROBS |
|-----|------------------|----|----|-----------|---------|
| 1   | chr1:7965399:G:A | A  | G  | 3.040e-04 | 401272  |

UKB_freq = UKB_freq %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
UKB_freq = UKB_freq %>% select(Start, MAF, A1, A2) %>% rename("UKB_MAF" = MAF) 
| Start   | UKB_MAF   | A1 | A2 |
|---------|-----------|----|----|
| 7965399 | 3.040e-04 | A  | G  |

mergename_freq3 = left_join(mergename_freq2, UKB_freq)
mergename_freq3 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq3)
# 964 83

---------------------------------------------
UKB_ALLPD = read.table("UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.assoc", header = T, sep = "")
| CHR | SNP              | BP      | A1 | F_A       | F_U      | A2 | CHISQ | P  | OR |
|-----|------------------|---------|----|-----------|----------|----|-------|----|----|
| 1   | chr1:7965399:G:A | 7965399 | A  | 0.0000000 | 0.000000 | G  | NA    | NA | NA |

UKB_ALLPD = UKB_ALLPD %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
UKB_ALLPD$End = UKB_ALLPD$Start
UKB_ALLPD = UKB_ALLPD %>% select(chr, Start, End, A1, A2, F_A, F_U, CHISQ, P, OR)
UKB_ALLPD = UKB_ALLPD %>% rename_with(~paste0(., "_UKB_Pheno_ALLPD_Freq"), F_A:OR)
| chr  | Start   | End     | A1 | A2 | F_A_UKB_Pheno_ALLPD_Freq | F_U_UKB_Pheno_ALLPD_Freq | CHISQ_UKB_Pheno_ALLPD_Freq | P_UKB_Pheno_ALLPD_Freq | OR_UKB_Pheno_ALLPD_Freq |
|------|---------|---------|----|----|--------------------------|--------------------------|----------------------------|------------------------|-------------------------|
| chr1 | 7965399 | 7965399 | A  | G  | 0.0000000                | 0.000000                 | NA                         | NA                     | NA                      |

mergename_freq3$End = as.character(mergename_freq3$End)
mergename_freq4 = left_join(mergename_freq3, UKB_ALLPD)

mergename_freq4 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq4)
# 964 88

---------------------------------------------
UKB_CASE_CTRL= read.table("UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC.assoc", sep = "", header = T)
# same layout as previous

UKB_CASE_CTRL = UKB_CASE_CTRL %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
UKB_CASE_CTRL$End = UKB_CASE_CTRL$Start
UKB_CASE_CTRL = UKB_CASE_CTRL %>% select(chr, Start, End, A1, A2, F_A, F_U, CHISQ, P, OR)
UKB_CASE_CTRL = UKB_CASE_CTRL %>% rename_with(~paste0(., "_UKB_Pheno_CASECTRL_Freq"), F_A:OR)

mergename_freq4$End = as.character(mergename_freq4$End)
mergename_freq5 = left_join(mergename_freq4, UKB_CASE_CTRL)

mergename_freq5 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq5)
# 964 93

---------------------------------------------
UKB_PARENT_CTRL = read.table("UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC.assoc", sep ="", header = T)
# layout same as previous

UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% separate(SNP, c("chr", "Start", "allele1", "allele2"), sep = ":")
UKB_PARENT_CTRL$End = UKB_PARENT_CTRL$Start
UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% select(chr, Start, End, A1, A2, F_A, F_U, CHISQ, P, OR)
UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% rename_with(~paste0(., "_UKB_Pheno_PARENTCTRL_Freq"), F_A:OR)

mergename_freq5$End = as.character(mergename_freq5$End)
mergename_freq6 = left_join(mergename_freq5, UKB_PARENT_CTRL)

mergename_freq6 %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 1 × 2
  Start        n
  <chr>    <int>
1 40322386     2

dim(mergename_freq6)
# 964 98

```

Re-arrange file
```
colnames(mergename_freq6)
mergename_freq6 = mergename_freq6 %>% select(Chr:assay.name, Gene.refGene, AAChange.refGene, pvalue:p.batch, AMP_MAF, N_INFORMATIVE_AMP_score:PVALUE_AMP_score, UKB_MAF, N_INFORMATIVE_UKB_ALL_score:OR_UKB_Pheno_PARENTCTRL_Freq) 

# clean AAChange.refGene column
mergename_freq6 = read.table("963_variants_score_AMP_UKB_wide.txt", header = T, sep = "\t")
mergename_freq6 %>% group_by(Gene.refGene) %>% tally() %>% arrange(desc(n)) %>% print(n=100)

# go through genes
mergename_freq6 %>% filter(Gene.refGene == "TRPM7") %>% select(Gene.refGene, AAChange.refGene)

# A tibble: 34 × 2
   Gene.refGene     n
   <chr>        <int>
 1 LRRK2          126 LRRK2:NM_198578
 2 GBA            116 PRKN:NM_004562 
 3 VPS13C         102 VPS13C:NM_020821
 4 POLG            97 POLG:NM_001126131, POLG:NM_002693
 5 PLA2G6          59 PLA2G6:NM_001004426, PLA2G6:NM_001349868
 6 PINK1           54 PINK1:NM_032409
 7 PRKN            54 PRKN:NM_004562
 8 ATP13A2         51 ATP13A2:NM_001141974
 9 EIF4G1          42 EIF4G1:NM_001194947
10 DNAJC13         41 DNAJC13:NM_001329126, DNAJC13:NM_015268
11 SYNJ1           41 SYNJ1:NM_003895
12 LRP10           25 LRP10:NM_014045
13 GIGYF2          24 GIGYF2:NM_001103148, GIGYF2:NM_001103146, GIGYF2:NM_001103147
14 FBXO7           22 FBXO7:NM_001033024, FBXO7:NM_012179
15 DNAJC6          21 DNAJC6:NM_001256864, DNAJC6:NM_001256865
16 HTRA2           10 HTRA2:NM_013247
17 PARK7           10 PARK7:NM_001123377, PARK7:NM_007262
18 VPS35           10 VPS35:NM_018206
19 TMEM230          9 TMEM230:NM_001330987, TMEM230:NM_001009923
20 SLC6A3           7 SLC6A3:NM_001044
21 SNCA             7 SNCA:NM_000345, SNCA:NM_001146054, SNCA:NM_001146055, SNCA:NM_001375286, SNCA:NM_001375288
22 SNCAIP           7 SNCAIP:NM_001308100
23 PINK1-AS         6 NONE (not actually of interest to us)
24 UCHL1            5 UCHL1:NM_004181
25 TNR              4 TNR:NM_003285
26 MAPT             3 MAPT:NM_001123067, MAPT:NM_001123066
27 RAB39B           2 RAB39B:NM_171998
28 SNCB             2 SNCB:NM_001318036, SNCB:NM_001318034, SNCB:NM_003085, SNCB:NM_001001502, SNCB:NM_001363140
29 TNK2             2 TNK2:NM_001308046, TNK2:NM_001382271, TNK2:NM_005781, TNK2:NM_001010938, TNK2:NM_001382272, TNK2:NM_001382273, TNK2:NM_001382274
30 FBXO7;FBXO7      1 FBXO7:NM_001033024, FBXO7:NM_012179
31 GLUD2            1 GLUD2:NM_012084
32 MRE11            1 MRE11:NM_001330347, MRE11:NM_005590, MRE11:NM_005591
33 NR4A2            1 NR4A2:NM_006186, NR4A2:NM_173173
34 TRPM7            1 TRPM7:NM_001301212, TRPM7:NM_017672

mergename_freq6 = mergename_freq6 %>% mutate(Gene.refGene = if_else(Gene.refGene == "FBXO7;FBXO7", "FBXO7", Gene.refGene))

mergename_freq6 = mergename_freq6 %>% mutate(AAChange.refGene = if_else(Gene.refGene == "LRRK2", "LRRK2:NM_198578",
                                                      if_else(Gene.refGene == "GBA", "PRKN:NM_004562",
                                                              if_else(Gene.refGene == "VPS13C", "VPS13C:NM_020821",
                                                                      if_else(Gene.refGene == "POLG", "POLG:NM_001126131, POLG:NM_002693",
                                                                              if_else(Gene.refGene == "PLA2G6", "PLA2G6:NM_001004426, PLA2G6:NM_001349868",
                                                                                      if_else(Gene.refGene == "PINK1", "PINK1:NM_032409",
                                                                                              if_else(Gene.refGene == "PRKN", "PRKN:NM_004562",
                                                                                                      if_else(Gene.refGene == "ATP13A2", "ATP13A2:NM_001141974",
                                                                                                              if_else(Gene.refGene == "EIF4G1", "EIF4G1:NM_001194947",
                                                                                                                      if_else(Gene.refGene == "DNAJC13", "DNAJC6:NM_001256864, DNAJC6:NM_001256865",
                                                                                                                              if_else(Gene.refGene == "SYNJ1", "SYNJ1:NM_003895",
                                                                                                                                      if_else(Gene.refGene == "LRP10", "LRP10:NM_014045",
                                                                                                                                              if_else(Gene.refGene == "GIGYF2", "GIGYF2:NM_001103148, GIGYF2:NM_001103146, GIGYF2:NM_001103147",
                                                                                                                                                      if_else(Gene.refGene == "FBXO7", "FBXO7:NM_001033024, FBXO7:NM_012179",
                                                                                                                                                              if_else(Gene.refGene == "DNAJC6", "DNAJC6:NM_001256864, DNAJC6:NM_001256865",
                                                                                                                                                                      if_else(Gene.refGene == "HTRA2", "HTRA2:NM_013247",
                                                                                                                                                                              if_else(Gene.refGene == "PARK7", "PARK7:NM_001123377, PARK7:NM_007262",
                                                                                                                                                                                      if_else(Gene.refGene == "VPS35", "VPS35:NM_018206",
                                                                                                                                                                                              if_else(Gene.refGene == "TMEM230", "TMEM230:NM_001330987, TMEM230:NM_001009923",
                                                                                                                                                                                                      if_else(Gene.refGene == "SLC6A3", "SLC6A3:NM_001044",
                                                                                                                                                                                                              if_else(Gene.refGene == "SNCA", "SNCA:NM_000345, SNCA:NM_001146054, SNCA:NM_001146055, SNCA:NM_001375286, SNCA:NM_001375288",
                                                                                                                                                                                                                      if_else(Gene.refGene == "SNCAIP", "SNCAIP:NM_001308100",
                                                                                                                                                                                                                              if_else(Gene.refGene == "UCHL1", "UCHL1:NM_004181",
                                                                                                                                                                                                                                      if_else(Gene.refGene == "TNR", "TNR:NM_003285",
                                                                                                                                                                                                                                              if_else(Gene.refGene == "MAPT", "MAPT:NM_001123067, MAPT:NM_001123066",
                                                                                                                                                                                                                                                      if_else(Gene.refGene == "RAB39B", "MAPT:NM_001123067, MAPT:NM_001123066",
                                                                                                                                                                                                                                                              if_else(Gene.refGene == "SNCB", "SNCB:NM_001318036, SNCB:NM_001318034, SNCB:NM_003085, SNCB:NM_001001502, SNCB:NM_001363140",
                                                                                                                                                                                                                                                                      if_else(Gene.refGene == "TNK2", "TNK2:NM_001308046, TNK2:NM_001382271, TNK2:NM_005781, TNK2:NM_001010938, TNK2:NM_001382272, TNK2:NM_001382273, TNK2:NM_001382274",
                                                                                                                                                                                                                                                                              if_else(Gene.refGene == "FBXO7", "FBXO7:NM_001033024, FBXO7:NM_012179",
                                                                                                                                                                                                                                                                                      if_else(Gene.refGene == "GLUD2", "GLUD2:NM_012084",
                                                                                                                                                                                                                                                                                              if_else(Gene.refGene == "MRE11", "MRE11:NM_001330347, MRE11:NM_005590, MRE11:NM_005591",
                                                                                                                                                                                                                                                                                                      if_else(Gene.refGene == "NR4A2", "NR4A2:NM_006186, NR4A2:NM_173173",
                                                                                                                                                                                                                                                                                                              if_else(Gene.refGene == "TRPM7", "TRPM7:NM_001301212, TRPM7:NM_017672", ""))))))))))))))))))))))))))))))))))
                                                                                                                                                                                                                                                                                                              
mergename_freq6 %>% group_by(Gene.refGene, AAChange.refGene) %>% tally() %>% print(n=100)

write.table(mergename_freq6, "963_variants_score_AMP_UKB_wide.txt", row.names = F, sep = "\t", quote = F)
```

DONE

