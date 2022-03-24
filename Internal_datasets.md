## Extract 963 23andme variants from internal data sets
- Authors: Vanessa Pitz, Cornelis Blauwendraat, Mary Makarious

## Objective
Extracting 963 variants (available in GRCH38 from 23andme data) from internal data sets and generate association stats on those variants only. Follow Mary's script for burden testing but instead of burden testing, do --score test Files can also be found here: https://drive.google.com/drive/u/1/folders/18kISsSBgz7iD9Gu9000tQupiY-SD1TOH Extracted variants from 23andme in GRCh37, AMP-PD are in GRCh38 rv Connect to biowulf and start node

```
sinteractive --mem=100g --cpus-per-task=20
```

# Workflow
0. Getting started
1. Generate 23andme variants file
2. Extract those variants from AMP-PDxNIH and UKB datasets
3. Association testing
4. Join all files together and write new file
5. Annotate file
6. Edit gene names by adding variant details


# 0. Getting started
### Working directory
```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
Create new variant file from old and new 23and me
```

# 1. Generate 23andme variants file
Combining two files we received
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

# right IDs
Ids = read.table("23andme_rsID_properID.txt", header = T, sep ="\t")
final = full_join(fulljoin_old_failed, Ids)

# old list
variants_1091 = read.table("List_1091variants_ID_only.txt", header = F, sep = "\t")

write.table(final, "results_collaborators_all_963.txt", row.names = F, quote = F, sep = "\t")

# create file with only IDs, so I can grep them from data sets
final = final %>% select(CHR.BP.REF.ALT)
names(final) = "V1"
head(final)

# check if these 23andme provided variants are all on our list
antijoin = anti_join(final, variants_1091)
dim(antijoin)
[1] 963   1

names(antijoin) = NULL
write.table(antijoin, "963_variants_in23andme_variants_IDs_only.txt", row.names = F, quote = F, sep = "\t")
q()
n
# All files now in GRCH38!
```

## UKB data contains indels

### PLINK file merged with no relateds

```
/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass.*
# these files are in plink2 format - convert

## .psam IDs
#IID	SEX
-000001	NA
-000002	NA
```

```
module load plink/2
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass --make-bed --out MERGED_UKB_first_pass

wc -l MERGED_UKB_first_pass.bim
# 16285684 MERGED_UKB_first_pass.bim
```

## Remove potential indels from UKB
### Write "no_indel" file

#### In bash
```
awk 'length($NF)==1 && length($(NF-1))==1' MERGED_UKB_first_pass.bim > MERGED_UKB_no_indels.txt
wc -l MERGED_UKB_no_indels.txt
# 14908659 MERGED_UKB_no_indels.txt

awk 'length($NF)>1 || length($(NF-1))>1' MERGED_UKB_first_pass.bim > UKB_indels.txt
wc -l UKB_indels.txt
# 1377025 UKB_indels.txt

cut -f 2 MERGED_UKB_no_indels.txt > UKB_noindels_tokeep.txt
```

#### In R
```
module load R
R
library(tidyverse)
library(data.table)

UKB = fread("MERGED_UKB_first_pass.bim")
dim(UKB)
# 16285684 
noindels = subset(UKB,nchar(V5)+nchar(V6)==2)
dim(noindels)
# 14908659
```

### Keep only non-indel variants
```
module load plink

plink --bfile MERGED_UKB_first_pass --extract UKB_noindels_tokeep.txt --make-bed --out MERGED_UKB_no_indels
# 14908659 variants and 200648 people pass filters and QC.
# Note: No phenotypes present.

# new working files: MERGED_UKB_no_indels
```


#### Double check if these are the correct variants
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
#828 1
q()
n
```


# 2. Extract 23andme variants from other datasets
## AMP-PD x LNG
File locations
```
AMP-PD
cd /data/CARD/PD/AMP_NIH/no_relateds/
# /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins*

#covariate file for burdens
/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt
```

### Grep variants
```
grep -w -f 963_variants_in23andme_variants_IDs_only.txt /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.bim > List_963variants_alleles_AMPPD.txt
wc -l List_963variants_alleles_AMPPD.txt
# 363 List_963variants_alleles_AMPPD.txt
```

### Double check in R
```
R
library(tidyverse)
library(data.table)

AMPPD = fread("/data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.bim")
togrep = fread("963_variants_in23andme_variants_IDs_only.txt"

AMMPD = AMPPD %>% dplyr::select(V2)
names(togrep) = "V2"

merge = merge(AMPPD, togrep)
dim(merge)
# 363 check
```

### Write new binary files
```
module load plink

plink --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins --extract List_963variants_alleles_AMPPD.txt --make-bed --out AMP_PD_963_23andme

# 363 variants and 7986 people pass filters and QC.
Among remaining phenotypes, 3376 are cases and 4610 are controls.
--make-bed to AMP_PD_963_23andme.bed + AMP_PD_963_23andme.bim +
AMP_PD_963_23andme.fam ... done.
AMP-PD only contains 363 out of 963 variants provided by 23andme
```

### Generate frequency file
```
module load plink
plink --bfile AMP_PD_963_23andme --freq --out AMP_PD_963_23andme

# 7986 people (4345 males, 3641 females) loaded from .fam.
# 7986 phenotype values loaded from .fam.
```

## UKB
### Grep variants
```
grep -w -f 963_variants_in23andme_variants_IDs_only.txt MERGED_UKB_no_indels.bim > List_963variants_alleles_UKB.txt
wc -l List_963variants_alleles_UKB.txt
# 718 List_963variants_alleles_UKB.txt
```

### Double check in R
```
R
library(tidyverse)
library(data.table)

UKB = fread("MERGED_UKB_no_indels.bim")
togrep = fread("963_variants_in23andme_variants_IDs_only.txt"

UKB = UKB %>% dplyr::select(V2)
names(togrep) = "V2"

merge = merge(UKB, togrep)
dim(merge)
# 718 check
```


### Write new binary files
```
module load plink/1.9.0-beta4.4

plink --bfile MERGED_UKB_no_indels --extract List_963variants_alleles_UKB.txt --make-bed --out UKB_963_23andme

# 718 variants and 200648 people pass filters and QC.
# Note: No phenotypes present.

# This file needs to be fixed because of an odd format where column 1 is all zeros (column 1 & 2 need to be same)
```

```
cut -d " " -f 2 UKB_963_23andme.fam > column2.txt
cut -d " " -f 2,3,4,5,6 UKB_963_23andme.fam > column23456.txt
scp UKB_963_23andme.fam UKB_963_23andme_ORIGINAL.fam
paste column2.txt column23456.txt > UKB_963_23andme.fam
```

### Generate frequency file
```
module load plink
plink --bfile UKB_963_23andme --freq --out UKB_963_23andme
```

# 3. Association testing
## Plink
### AMP-PD x LNG
```
plink --bfile AMP_PD_963_23andme --assoc --pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --pheno-name PD_PHENO --out AMP_PD_963_23andme_pheno

# 363 variants and 7986 people pass filters and QC.
# Among remaining phenotypes, 3376 are cases and 4610 are controls.
# Writing C/C --assoc report to AMP_PD_963_23andme_pheno.assoc ...
done.
```

### UKB
#### ALL_PD
```
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC

# 718 variants and 200648 people pass filters and QC.
# Among remaining phenotypes, 7806 are cases and 38051 are controls.  (154791
phenotypes are missing.)
# Writing C/C --assoc report to
# UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.assoc ...
done.

```

#### PD_CASE_CONTROL
```
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC

# 718 variants and 200648 people pass filters and QC.
# Among remaining phenotypes, 1105 are cases and 5643 are controls.  (193900
phenotypes are missing.)
# Writing C/C --assoc report to
# UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC.assoc ...
done.

```

#### PD_PARENT_CONTROL_2021_with_PC
```
plink --bfile UKB_963_23andme --assoc --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --pheno-name PHENO --allow-no-sex --out UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC

# 718 variants and 200648 people pass filters and QC.
# Among remaining phenotypes, 6033 are cases and 28945 are controls.  (165670
phenotypes are missing.)
# Writing C/C --assoc report to
# UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC.assoc ...
done.

```
## Rvtest
### AMP-PD x LNG
#### Convert binary files to vcf for input
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

#### Run Rvtest
```
module load rvtests

rvtest --inVcf AMP_PD_963_23andme.vcf.gz --pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --pheno-name PD_PHENO --out AMP_PD_963_23andme_withcovars_score --single wald,score --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5

# new files
AMP_PD_963_23andme_withcovars_score.log
AMP_PD_963_23andme_withcovars_score.SingleScore.assoc
AMP_PD_963_23andme_withcovars_score.SingleWald.assoc

# Pheno file
AMP_pheno %>% group_by(PD_PHENO) %>% tally()
# A tibble: 2 × 2
  PD_PHENO     n
     <int> <int>
1        1  4610 
2        2  3376

(1 = CONTROL, 2 = CASE)
```

### UKB
#### Convert binary files to vcf for input
```
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


#### ALL_PD
```
# Phenotype file looks like this
# using UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
# # A tibble: 4 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL    38051
2 parent      6033
3 PD          1105
4 sibling      668

# UKB_all %>% group_by(SEX) %>% tally()
# A tibble: 2 × 2
    SEX     n
  <int> <int>
1     1 22040
2     2 23817


# now run rvtest
module load rvtests

rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files generated
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleWald.assoc
```

#### PD_CASE_CONTROL
```
# Phenotype file looks like this
# using UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
# # A tibble: 2 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL     5643
2 PD          1105

# run rvtest

rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleWald.assoc
```

#### PD_PARENT_CONTROL_2021_with_PC
```
# Phenotype file looks like this
# using UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
# # A tibble: 2 × 2
  PHENO_NAME     n
  <chr>      <int>
1 CONTROL    28945
2 parent      6033

# run rvtest
rvtest --inVcf UKB_963_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --pheno-name PHENO --out UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC --single wald,score --covar /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt --covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5

## new files 
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.log
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleScore.assoc
UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleWald.assoc
```

# 4. Generate final file
### Read all files into R
Change layout of files to match

#### 23andme variants file
```
R
library(dplyr)
library(tidyr)


#####
andme = read.table("results_collaborators_all_963.txt", sep = "\t", header = T)
dim(andme)

# layout
| assay.name  | scaffold | position  | alleles |CHR.BP.REF.ALT|
|-------------|----------|-----------|---------|--------------|
| rs536425950 | chr1     | 155235024 | C/T     | chr1:155235024:T:C |

# change
andme = andme %>% separate(CHR.BP.REF.ALT, c("Chr", "Start", "REF", "ALT"), sep = ":")
andme$End = andme$Start
andme = andme %>% select(Chr, Start, End, REF, ALT, assay.name, pvalue:p.batch)
dim(andme)
# 963  35
```


#### 23andme variants found in AMP-PD x LNG 
```
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
# 349  14

# layout we need for annotation
| Chr  | Start     | End       | Ref | Alt |
|------|-----------|-----------|-----|-----|
| chr1 | 155235006 | 155235006 | A   | G   |

```


#### 23andme variants found in UKB
All PD
```
UKB_ALL_score = read.table("UKB_963_23andme_withcovars_score_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.SingleScore.assoc",  sep = "\t", header = T)

# change
UKB_ALL_score = UKB_ALL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_ALL_score$End = UKB_ALL_score$Start
UKB_ALL_score$Chr = paste0('chr', UKB_ALL_score$Chr)
UKB_ALL_score = UKB_ALL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_ALL_score = UKB_ALL_score %>% rename_with(~paste0(., "_UKB_ALL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_ALL_score)
# 718  14
```

Case control
```
UKB_CASE_CTRL_score = read.table("UKB_963_23andme_withcovars_score_PD_CASE_CONTROL_2021_with_PC.SingleScore.assoc", sep = "\t", header = T)

# change
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_CASE_CTRL_score$End = UKB_CASE_CTRL_score$Start
UKB_CASE_CTRL_score$Chr = paste0('chr', UKB_CASE_CTRL_score$Chr)
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_CASE_CTRL_score = UKB_CASE_CTRL_score %>% rename_with(~paste0(., "_UKB_CASE_CONTROL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_CASE_CTRL_score)
# 718  14
```

Parent control
```
UKB_PARENT_CTRL_score = read.table("UKB_963_23andme_withcovars_score_PD_PARENT_CONTROL_2021_with_PC.SingleScore.assoc", sep = "\t", header = T)

# change
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% rename("Chr" = CHROM, "Start" = POS)
UKB_PARENT_CTRL_score$End = UKB_PARENT_CTRL_score$Start
UKB_PARENT_CTRL_score$Chr = paste0('chr', UKB_PARENT_CTRL_score$Chr)
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% select(Chr, Start, End, REF, ALT, N_INFORMATIVE:PVALUE)
UKB_PARENT_CTRL_score = UKB_PARENT_CTRL_score %>% rename_with(~paste0(., "_UKB_PARENT_CONTROL_score"), N_INFORMATIVE:PVALUE)
dim(UKB_PARENT_CTRL_score)
# 718  14
```

#### Start merging
```
andme$Start =as.integer(andme$Start)
andme$End =as.integer(andme$End)

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
# if all FALSE, then this is good

write.table(merge4, "963variants_23andme_AMPPD_UKB_cov_merged_score.txt", row.names = F, sep = "\t", quote= F)
```


#### Write variant file for annotation
```
R 
library(tidyverse)
to_annotate = merge4 %>% select(Chr, Start, End, REF, ALT)
names(to_annotate) = NULL
write.table(to_annotate, "963variants_to_annotate.txt", sep = "\t", row.names =F, quote= F)
```

# 5. Annotation
### Run ANNOVAR
Data is in genome build hg38
```
module load annovar

#gene build hg38

table_annovar.pl 963variants_to_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150 \
-operation g,f -outfile 963variants_23_AMP_UKB -nastring .

## new files
963variants_23_AMP_UKB.hg38_multianno.txt
```

### Join gene names onto final 
#### Edit file
```
R
library(dplyr)
library(tidyr)

variants = read.table("963variants_23andme_AMPPD_UKB_cov_merged_score.txt", header = T, sep = "\t")
variants = variants %>% tidyr::unite("CHR.START.REF.ALT", c("Chr", "Start", "REF", "ALT"), sep = ":")

Genenames = read.table("963variants_23_AMP_UKB.hg38_multianno.txt", header = T, sep = "\t")
Genenames = Genenames %>% tidyr::unite("CHR.START.REF.ALT", c("Chr", "Start", "Ref", "Alt"), sep = ":")

Genenames = Genenames %>% select(CHR.START.REF.ALT, Gene.refGene, AAChange.refGene)

#merge Gene name onto list
mergename = left_join(variants, Genenames)
dim(mergename)
# 963 70
```

#### Add frequency files
```
# MAF for AMP
AMP_freq = read.table("AMP_PD_963_23andme.frq", header = T, sep = "")
| Chr | SNP              | A1 | A2 | MAF       | NCHROBS |
|-----|------------------|----|----|-----------|---------|
| 1   | chr1:7965399:G:A | A  | G  | 0.0001878 | 15972   |

## change layout
AMP_freq = AMP_freq %>% rename("CHR.START.REF.ALT"=SNP)
AMP_freq = AMP_freq %>% select(CHR.START.REF.ALT, MAF) %>% rename("AMP_MAF" = MAF)

mergename$Start = as.character(mergename$Start)
mergename_freq1 = left_join(mergename, AMP_freq)
mergename_freq1 %>% group_by(CHR.START.REF.ALT) %>% tally() %>% filter(n>1)

dim(mergename_freq1)
# 963 71
```

```
# MAF for phenotypes
AMP_phenoMAF = read.table("AMP_PD_963_23andme_pheno.assoc", header = T, sep = "")

| CHR | SNP              | BP      | A1 | F_A       | F_U       | A2 | CHISQ     | P      | OR     |
|-----|------------------|---------|----|-----------|-----------|----|-----------|--------|--------|
| 1   | chr1:7965399:G:A | 7965399 | A  | 0.0001481 | 0.0002169 | G  | 9.829e-02 | 0.7539 | 0.6827 |

AMP_phenoMAF = AMP_phenoMAF %>% rename("CHR.START.REF.ALT"=SNP)
AMP_phenoMAF = AMP_phenoMAF %>% select(CHR.START.REF.ALT, F_A, F_U, CHISQ, P, OR)
AMP_phenoMAF = AMP_phenoMAF %>% rename_with(~paste0(., "_AMP_PhenoFreq"), F_A:OR)

| chr  | Start   | End     | A1 | A2 | F_A_AMP_PhenoFreq | F_U_AMP_PhenoFreq | CHISQ_AMP_PhenoFreq | P_AMP_PhenoFreq | OR_AMP_PhenoFreq |
|------|---------|---------|----|----|-------------------|-------------------|---------------------|-----------------|------------------|
| chr1 | 7965399 | 7965399 | A  | G  | 0.0001481         | 0.0002169         | 9.829e-02           | 0.7539          | 0.6827           |
```

Merge
```
mergename_freq1$End = as.character(mergename_freq1$End)
mergename_freq2 = left_join(mergename_freq1, AMP_phenoMAF)

mergename_freq2 %>% group_by(CHR.START.REF.ALT) %>% tally() %>% filter(n>1)

dim(mergename_freq2)
# 963 76
```

Add UKB frequency (phenotype and MAF) files
```
UKB_freq = read.table("UKB_963_23andme.frq", header = T, sep = "")

| CHR | SNP              | A1 | A2 | MAF       | NCHROBS |
|-----|------------------|----|----|-----------|---------|
| 1   | chr1:7965399:G:A | A  | G  | 3.040e-04 | 401272  |

UKB_freq = UKB_freq %>% rename("CHR.START.REF.ALT"=SNP)

UKB_freq = UKB_freq %>% select(CHR.START.REF.ALT, MAF) %>% rename("UKB_MAF" = MAF) 

| Start   | UKB_MAF   | A1 | A2 |
|---------|-----------|----|----|
| 7965399 | 3.040e-04 | A  | G  |

mergename_freq3 = left_join(mergename_freq2, UKB_freq)
mergename_freq3 %>% group_by(CHR.START.REF.ALT) %>% tally() %>% filter(n>1)

dim(mergename_freq3)
# 963 77
```

Add UKB_ALL_PD.assoc
```
UKB_ALLPD = read.table("UKB_963_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.assoc", header = T, sep = "")

| CHR | SNP              | BP      | A1 | F_A       | F_U      | A2 | CHISQ | P  | OR |
|-----|------------------|---------|----|-----------|----------|----|-------|----|----|
| 1   | chr1:7965399:G:A | 7965399 | A  | 0.0000000 | 0.000000 | G  | NA    | NA | NA |

UKB_ALLPD = UKB_ALLPD %>% rename("CHR.START.REF.ALT"=SNP)

UKB_ALLPD = UKB_ALLPD %>% select(CHR.START.REF.ALT, F_A, F_U, CHISQ, P, OR)
UKB_ALLPD = UKB_ALLPD %>% rename_with(~paste0(., "_UKB_Pheno_ALLPD_Freq"), F_A:OR)

| chr  | Start   | End     | A1 | A2 | F_A_UKB_Pheno_ALLPD_Freq | F_U_UKB_Pheno_ALLPD_Freq | CHISQ_UKB_Pheno_ALLPD_Freq | P_UKB_Pheno_ALLPD_Freq | OR_UKB_Pheno_ALLPD_Freq |
|------|---------|---------|----|----|--------------------------|--------------------------|----------------------------|------------------------|-------------------------|
| chr1 | 7965399 | 7965399 | A  | G  | 0.0000000                | 0.000000                 | NA                         | NA                     | NA                      |
```

Merge
```
mergename_freq3$End = as.character(mergename_freq3$End)
mergename_freq4 = left_join(mergename_freq3, UKB_ALLPD)

dim(mergename_freq4)
# 964 82
```

Add UKB_CASE_CTRL.assoc
```
UKB_CASE_CTRL= read.table("UKB_963_23andme_PD_CASE_CONTROL_2021_with_PC.assoc", sep = "", header = T)
# same layout as previous

UKB_CASE_CTRL = UKB_CASE_CTRL %>% rename("CHR.START.REF.ALT"=SNP)
UKB_CASE_CTRL = UKB_CASE_CTRL %>% select(CHR.START.REF.ALT, F_A, F_U, CHISQ, P, OR)
UKB_CASE_CTRL = UKB_CASE_CTRL %>% rename_with(~paste0(., "_UKB_Pheno_CASECTRL_Freq"), F_A:OR)
```

```
mergename_freq4$End = as.character(mergename_freq4$End)
mergename_freq5 = left_join(mergename_freq4, UKB_CASE_CTRL)

dim(mergename_freq5)
# 963 87
```

Add UKB_PARENT_CTRL.assoc
```
UKB_PARENT_CTRL = read.table("UKB_963_23andme_PD_PARENT_CONTROL_2021_with_PC.assoc", sep ="", header = T)
# layout same as previous

UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% rename("CHR.START.REF.ALT"=SNP)
UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% select(CHR.START.REF.ALT, F_A, F_U, CHISQ, P, OR)
UKB_PARENT_CTRL = UKB_PARENT_CTRL %>% rename_with(~paste0(., "_UKB_Pheno_PARENTCTRL_Freq"), F_A:OR)
```

Merge
```
mergename_freq5$End = as.character(mergename_freq5$End)
mergename_freq6 = left_join(mergename_freq5, UKB_PARENT_CTRL)

dim(mergename_freq6)
# 963 92
```

# 6. Edit gene names by adding variant details
```
colnames(mergename_freq6)
mergename_freq6 = mergename_freq6 %>% dplyr::select(CHR.START.REF.ALT:assay.name, Gene.refGene, AAChange.refGene, pvalue:p.batch, AMP_MAF, N_INFORMATIVE_AMP_score:PVALUE_AMP_score, UKB_MAF, N_INFORMATIVE_UKB_ALL_score:OR_UKB_Pheno_PARENTCTRL_Freq) 

# write file
write.table(mergename_freq6, "963_variants_score_AMP_UKB_AAChangenotclean.txt", quote = F, sep = "\t", row.names = F)

# clean AAChange.refGene column
# mergename_freq6 = read.table("963_variants_score_AMP_UKB_AAChangenotclean.txt", header = T, sep = "\t")
mergename_freq6 %>% group_by(Gene.refGene) %>% tally() %>% arrange(desc(n)) %>% print(n=100)

# go through genes
mergename_freq6 %>% filter(Gene.refGene == "TRPM7") %>% select(Gene.refGene, AAChange.refGene)

# A tibble: 34 × 2
   Gene.refGene     n
   <chr>        <int>
 1 LRRK2          125 LRRK2:NM_198578
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

Move on to meta-analysis.
