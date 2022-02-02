## Extract 1091 variants from internal data sets
- Vanessa

---
### Quick Description: 
Extracting 1091 variants (derived from 23andme data) from internal data sets and generate association stats on those variants only.
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

Now take those files and merge them with the latest 23andme files
```
R
library(dplyr)
library(tidyr)


#####
andme = read.table("results_collaborators_all_963.txt", sep = ",", header = T)

# layout
| assay.name  | scaffold | position  | alleles |
|-------------|----------|-----------|---------|
| rs536425950 | chr1     | 155235024 | C/T     |

# change
andme = andme %>% rename("Chr" = scaffold, "Start" = position)
andme$End = andme$Start
andme = andme %>% separate(alleles, c("REF", "ALT"), sep = "/")
andme = andme %>% select(Chr, Start, End, REF, ALT, assay.name, pvalue:p.batch)


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


#####
AMP_wald = read.table("AMP_PD_963_23andme_withcovars_score.SingleWald.assoc", sep = "\t", header = T)
# layout 
| CHROM | POS     | REF | ALT | N_INFORMATIVE |
|-------|---------|-----|-----|---------------|
| 1     | 7965399 | G   | A   | 7974          |

AMP_wald = AMP_wald %>% rename("Chr" = CHROM, "Start" = POS)
AMP_wald$End = AMP_wald$Start
AMP_wald$Chr = paste0('chr', AMP_wald$Chr)
AMP_wald = AMP_wald %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:Pvalue)
AMP_wald = AMP_wald %>% rename_with(~paste0(., "_AMP_wald"), N_INFORMATIVE:Pvalue)
AMP_wald = AMP_wald %>% rename("Test" = Test_AMP_wald)

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


#####
UKB_ALL_wald = read.table("UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleWald.assoc",  sep = "\t", header = T)

# change
UKB_ALL_wald = UKB_ALL_wald %>% rename("Chr" = CHROM, "Start" = POS)
UKB_ALL_wald$End = UKB_ALL_wald$Start
UKB_ALL_wald$Chr = paste0('chr', UKB_ALL_wald$Chr)
UKB_ALL_wald = UKB_ALL_wald %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:Pvalue)
UKB_ALL_wald = UKB_ALL_wald %>% rename_with(~paste0(., "_UKB_ALL_wald"), N_INFORMATIVE:Pvalue)
UKB_ALL_wald = UKB_ALL_wald %>% rename("Test" = Test_UKB_ALL_wald)
```

Start merging
```
merge1 = left_join(andme, AMP_score)
dim(merge1)
# 963 44


merge2 = left_join(merge1, AMP_wald)
dim(merge2)
# 2111 49


merge3 = left_join(merge2, UKB_ALL_score)
merge4 = left_join(merge3, UKB_ALL_wald)

```




