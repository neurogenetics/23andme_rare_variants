## Extract 1091 variants from internal data sets
- Vanessa

---
### Quick Description: 
Extracting 1091 variants (derived from 23andme data) from internal data sets and generate association stats on those variants only.
Follow Mary's script for burden testing but instead of burden testing, do --score test

Connect to biowulf and start node
```
sinteractive --mem=100g --cpus-per-task=20
# new folder
/data/CARD/projects/23andme_annotation/variants_in_internal_datasets

# variant file sits here
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
/data/CARD/projects/23andme_annotation/variants_in_internal_datasets/List_1091variants_ID2.txt 
```

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
grep -w -f List_1091variants_ID2.txt /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.bim > List_1091variants_alleles_AMPPD.txt

wc -l List_1091variants_alleles_AMPPD.txt
# 19 List_1091variants_alleles_AMPPD.txt
```

Double check if these are the correct variants
```
R
library(dplyr)
library(tidyr)
variants_new = read.table("List_1091variants_alleles_AMPPD.txt", header =F, sep = "\t")
variants = read.table("List_1091variants_ID2.txt",  header =F, sep = "\t")

variants_new = variants_new %>% separate(V2, c("chr", "position", "allele1", "allele2"), sep = ":")
Check = variants_new %>% select(chr, position)
Check = Check %>% unite("V1",c(chr:position), sep = ":")
join = left_join(Check, variants)
dim(join)
# 19 1
q()
n
```

Write new binary files - only 19 variants
```
module load plink

plink --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins --extract List_1091variants_alleles_AMPPD.txt --make-bed --out AMP_PD_1091_23andme


# 386449 MB RAM detected; reserving 193224 MB for main workspace.
# 47899791 variants loaded from .bim file.
# 7986 people (4345 males, 3641 females) loaded from .fam.
# 7986 phenotype values loaded from .fam.
# --extract: 19 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 7986 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# 19 variants and 7986 people pass filters and QC.
# Among remaining phenotypes, 3376 are cases and 4610 are controls.
# --make-bed to AMP_PD_1091_23andme.bed + AMP_PD_1091_23andme.bim +
# AMP_PD_1091_23andme.fam ... done.
```
AMP-PD only contains 19 out of 1091 variants found in 23andme

Run association test in rvtest
```
# convert plink files to vcf for rvtest
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile AMP_PD_1091_23andme \
--export vcf bgz id-paste=iid --out AMP_PD_1091_23andme${OUTNAME} --mac 1

tabix -p vcf  AMP_PD_1091_23andme${OUTNAME}.vcf.gz

# new files
AMP_PD_1091_23andme.vcf.gz
AMP_PD_1091_23andme.vcf.gz.tbi
```

Run rvtests
```
module load rvtests
rvtest --inVcf AMP_PD_1091_23andme.vcf --pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --out AMP_PD_1091_23andme_score --score

# Unparsed arguments:  --score
```

### UKB
```
# PLINK file merged with no relateds
/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass.*
# these files are in plink2 format - convert

module load plink/2
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass --make-bed --out MERGED_UKB_first_pass

grep -w -f List_1091variants_ID2.txt MERGED_UKB_first_pass.bim > List_1091variants_alleles_UKB.txt
wc -l List_1091variants_alleles_UKB.txt
# 19 List_1091variants_alleles_UKB.txt
```

Double check if these are the correct variants
```
R
library(dplyr)
library(tidyr)
variants_new = read.table("List_1091variants_alleles_UKB.txt", header =F, sep = "\t")
variants = read.table("List_1091variants_ID2.txt",  header =F, sep = "\t")

variants_new = variants_new %>% separate(V2, c("chr", "position", "allele1", "allele2"), sep = ":")
Check = variants_new %>% dplyr::select(chr, position)
Check = Check %>% unite("V1",c(chr:position), sep = ":")
join = left_join(Check, variants)
dim(join)
# 19 1
q()
n
```

Write new binary files - only 19 variants

```
module load plink/1.9.0-beta4.4

plink --bfile MERGED_UKB_first_pass --extract List_1091variants_alleles_UKB.txt --make-bed --out UKB_1091_23andme

# 386449 MB RAM detected; reserving 193224 MB for main workspace.
# 16285684 variants loaded from .bim file.
# 200648 people (0 males, 0 females, 200648 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to UKB_1091_23andme.nosex .
# --extract: 19 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 200648 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.997511.
# 19 variants and 200648 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to UKB_1091_23andme.bed + UKB_1091_23andme.bim +
# UKB_1091_23andme.fam ... done.
```

Run association test in rvtest
```
# convert binary files to vcf for input in rvtest
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile UKB_1091_23andme \
--export vcf bgz id-paste=iid --out UKB_1091_23andme${OUTNAME} --mac 1

tabix -p vcf  UKB_1091_23andme${OUTNAME}.vcf.gz

# new files written: 
UKB_1091_23andme.vcf.gz
UKB_1091_23andme.vcf.gz.tbi
```

Run rvtests
```
module load rvtests

rvtest --inVcf UKB_1091_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --out UKB_1091_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC --score

####
[pitzv2@cn0896 variants_in_internal_datasets]$ rvtest --inVcf UKB_1091_23andme.vcf.gz --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt --out UKB_1091_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC --score
Thank you for using rvtests (version: 20190205, git: c86e589efef15382603300dc7f4c3394c82d69b8)
  For documentations, refer to http://zhanxw.github.io/rvtests/
  For questions and comments, plase send to Xiaowei Zhan <zhanxw@umich.edu>
  For bugs and feature requests, please submit at: https://github.com/zhanxw/rvtests/issues

The following parameters are available.  Ones with "[]" are in effect:

Available Options
      Basic Input/Output: --inVcf [UKB_1091_23andme.vcf.gz], --inBgen []
                          --inBgenSample [], --inKgg []
                          --out [UKB_1091_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC]
                          --outputRaw
       Specify Covariate: --covar [], --covar-name [], --sex
       Specify Phenotype:
                          --pheno [/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt]
                          --inverseNormal, --useResidualAsPhenotype, --mpheno []
                          --pheno-name [], --qtl, --multiplePheno []
        Specify Genotype: --dosage [], --multipleAllele
    Chromosome X Options: --xLabel [], --xParRegion []
           People Filter: --peopleIncludeID [], --peopleIncludeFile []
                          --peopleExcludeID [], --peopleExcludeFile []
             Site Filter: --rangeList [], --rangeFile [], --siteFile []
                          --siteDepthMin [], --siteDepthMax [], --siteMACMin []
                          --annoType []
         Genotype Filter: --indvDepthMin [], --indvDepthMax [], --indvQualMin []
       Association Model: --single [], --burden [], --vt [], --kernel []
                          --meta []
     Family-based Models: --kinship [], --xHemiKinship [], --kinshipEigen []
                          --xHemiKinshipEigen [], --boltPlink []
                          --boltPlinkNoCheck
          Grouping Unit : --geneFile [], --gene [], --setList [], --setFile []
                          --set []
        Frequency Cutoff: --freqUpper [], --freqLower []
            Missing Data: --impute [], --imputePheno, --imputeCov
    Conditional Analysis: --condition []
     Auxiliary Functions: --noweb, --hide-covar, --numThread [], --outputID
                          --help

Effective Options
    --inVcf UKB_1091_23andme.vcf.gz
    --out UKB_1091_23andme_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC
    --pheno /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt

# Covariate files
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt

```





