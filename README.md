## Re-do variant selection 
 - **Project:** Assessment of Parkinson's disease pathogenic genes and risk variants
 - **Author(s):** Cornelis, Mike, Andy, and Vanessa
 - **Date Last Updated:** July 2022

### Quick Description: 
Multiple genes and rare variants are associated with Parkinson’s disease and have been shown to be either a strong risk factor or to cause disease. However, the majority of genes are under heavy debate about their actual pathogenicity and this analysis aims to assess what genes and rare variants are actually associated with Parkinson’s disease.

### Motivation/Goals:
Rare variant analyses in known PD genes
- Create a list of variants from 23andMe arrays that are reported to be associated in the literature with Parkinson’s disease (or related disease) => note this can be done by us

- Perform rare variant case-control association testing per sub-population (optional, but preferred) of each variant of interest. Currently for rare variant testing we use “Score test” but we are open for suggestions on 23andMe preferred tests.

- Investigate these findings and meta-analyze in other datasets such as IPDGC (https://pdgenetics.org/), GP2 (https://www.parkinsonsroadmap.org/gp2/) and UKbiobank data.

- Optional, perform burden testing in genes using variants of interest using standard burden test algorithms like SKAT, CMC

- Report these findings in collaboration with the 23andMe Team in a scientific journal.

## Structure of README:

### [1. Understanding data and subset data](#1-understanding-the-data)
This section goes through:
- understanding the data and where is the data located

### [2. Data cleaning and annotation](#2-data-cleaning-and-annotation)
How to clean the data and and annotate

### [3. Extract gene-specific information](#3-gene-specific-information)
Extracting gene-specific information especially with focus on genes identified in PD review paper Blauwendraat et al 2019 (PMI: 31521533) 

### [4. Generate final list of variants](#4-final-list-of-variants)
Write a final file with variants of interest from this data set


## 1. Understanding the data
Data can be located in /data/CARD/projects/23andme_annotation. Note: the data is in genome build hg19.
All SNPs are in the file "all_snp_info.txt".

```
#get node
sinteractive --mem=100g --cpus-per-task=20

cd /data/CARD/projects/23andme_annotation
wc -l all_snp_info_txt

# 64,523,260 all_snp_info.txt
```

### Create new folder 
```
mkdir mkdir Redo_variants_July2022
```

Make short file
```
cut -f 4,5,6,7,10 all_snp_info.txt > /data/CARD/projects/23andme_annotation/Redo_variants_July2022/short_all_snp_info.txt

cd ./Redo_variants_July2022
```

```
| assay.name | scaffold | position | alleles | gene.context |
|------------|----------|----------|---------|--------------|
```


### General idea:
- ignore indels since no alleles present :(
- remove multi-allelics since its always tricky with them :(
- annotate with annovar using 23andme alleles
- extract all coding regions
- get correct hg19 REF alleles
- flip alleles to make them right
- re-annotate with annovar with correct alleles
- extract genes of interest


## 2. Data cleaning and annotation

D/I for alleles indicates that they are missing here.
```
#only writing columns 'scaffold' and 'position'
cut -f 2,3 short_all_snp_info.txt > temp1.txt 

#only writing column 'position'
cut -f 3 short_all_snp_info.txt > temp2.txt

#only writing column 'alleles' and tab-separate them. e.g. A/G -> A  G
cut -f 4 short_all_snp_info.txt | sed 's/\//\t/g' > temp3.txt

#only writing column 'assay.name' and 'gene.context'
cut -f 1,5 short_all_snp_info.txt > temp4.txt

#combining the temp files, to end up with a new format
paste temp1.txt temp2.txt temp3.txt temp4.txt > input_annovar_first_pass.txt

| scaffold | position | position | alleles | assay.name | gene.context |
|----------|----------|----------|---------|------------|--------------|
```

```
wc -l input_annovar_first_pass.txt
# 64,523,260 input_annovar_first_pass.txt
```

First annotation with ANNOVAR
```
module load annovar

table_annovar.pl input_annovar_first_pass.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150 \
-operation g,f -outfile first_pass_23andme -nastring .

# This step creates the file "first_pass_23andme.hg19_multianno.txt", which looks like this:
| Chr | Start | End | Ref | Alt | Func.refGene | Gene.refGene | GeneDetail.refGene | ExonicFunc.refGene | AAChange.refGene | avsnp150 |
|-----|-------|-----|-----|-----|--------------|--------------|--------------------|--------------------|------------------|----------|

refGene from ANNOVAR adds these columns: Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, while avsnp150 adds rsnumbers using dbSNP annotations. 

(just ran this on the side in a swarm)
swarm -f /data/CARD/projects/23andme_annotation/Redo_variants_July2022/annotate_1.swarm --verbose 1 --module annovar -g 100 -t auto 
```

```
wc -l first_pass_23andme.hg19_multianno.txt 
# 64,523,261 first_pass_23andme.hg19_multianno.txt
# still same number as before - check
```

Check how many unique rows 
```
sort first_pass_23andme.hg19_multianno.txt | uniq | wc -l
# 64,422,687 (with header)
```

Double-check in R
```
R
library(data.table)
data = fread("first_pass_23andme.hg19_multianno.txt")
data %>% distinct() %>% tally()
# 64,422,686
q()
n
```

Create new file with only unique rows and variant identifier, combining columns Chr and Start, separating them with :
```
cut -f 1,2 first_pass_23andme.hg19_multianno.txt | sed 's/\t/:/g' | sort | uniq -d > multiallelics.txt

wc -l multiallelics.txt 
# File 671,848 multiallelics.txt (no header) - 671,848 multiallelic variants to be removed
```

### Check D/I variants

4,124,008 variants have alleles labelled as D/I but are they useful?
Write them into a new file by replacing D with A and I with G and then annotating again with ANNOVAR:

```
awk '$4 == "D"' input_annovar_first_pass.txt > input_annovar_deletions.txt
wc -l input_annovar_deletions.txt 
# 4,124,008 input_annovar_deletions.txt (no header)

#replace the alleles
sed -i -e 's/D/A/g' input_annovar_deletions.txt
sed -i -e 's/I/G/g' input_annovar_deletions.txt
```

Now annotate those to see if they're interesting
```
table_annovar.pl input_annovar_deletions.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150 \
-operation g,f -outfile input_annovar_deletions_annovar_AG_converted -nastring .

#new file written input_annovar_deletions_annovar_AG_converted.hg19_multianno.txt
```

Annotation looks like these variants could be useful but without alleles, the data can not be reliably analysed. So further analysis is based on alleles A-C-T-G and ignores D/I.

Extract 907,876 exonic variants and write into a new file
```
grep "exonic" first_pass_23andme.hg19_multianno.txt > exonic_only.txt
wc -l exonic_only.txt
# 907,876 exonic_only.txt (no header)
```

Extract 11,667 splicing variants and write into a new file
```
grep "splicing" first_pass_23andme.hg19_multianno.txt > splicing_only.txt
wc -l splicing_only.txt 
# 11,667 splicing_only.txt (no header)
```

Now merge exonic_only.txt and splicing_only.txt
```
cat exonic_only.txt splicing_only.txt | sort -u > exonic_splicing_only.txt
wc -l exonic_splicing_only.txt
# 918,876 exonic_splicing_only.txt
```

And add variant name
```
#variant name style chr10:100003915
cut -f 1,2 exonic_splicing_only.txt | sed 's/\t/:/g' > var_name.txt

paste exonic_splicing_only.txt var_name.txt > exonic_splicing_only_v2.txt

wc -l exonic_splicing_only_v2.txt
# 918,876 exonic_splicing_only_v2.txt
```

Fixing space issue
```
sed -i 's/ /_/g' exonic_splicing_only_v2.txt
# this removes the spaces in column 9, e.g. nonsynonymous SNV -> nonsynonymous_SNV
```

### Remove multi-allelic variants
```
module load R/3.6.0 
R
require(data.table)
library(dplyr)
data <- fread("exonic_splicing_only_v2.txt",header=F)
dim(data)
# 918876     12

dup <- fread("multiallelics.txt",header=F)
names(dup)[1] <- "V12"
dim(dup)
# 671848      1


# anti-join here
data_no_dup <- anti_join(data, dup, by="V12")
dim(data_no_dup)
# 903097     12

double_check <- merge(data, dup, by="V12")
dim(double_check)
# 15779    12

write.table(data_no_dup, file="exonic_splicing_only_no_dup.txt", quote=FALSE,row.names=F,sep="\t")
q()
n
```

```
wc -l exonic_splicing_only_no_dup.txt
# 903,098 exonic_splicing_only_no_dup.txt
```

Duplicate variants that were present in dup *and* in data were removed. Since the multiallelics.txt and exonix_splicing_only_v2.txt files are both based on first_pass_23andme.hg19_multianno.txt, the number of rows left in the file makes sense.
A total of 903,097 variants are included in this file.

Get correct reference allele
```
cut -f 1,2 exonic_splicing_only_no_dup.txt | grep -v "V1" > prep_for_ref_allele.txt

## Add headers to file
echo "chr" > header.txt
# "chr" column header
echo "start" > header1.txt
# "start" column header
paste header.txt header1.txt > header_final.txt
# "chr" "start" column
cat header_final.txt prep_for_ref_allele.txt > prep_for_ref_allele_input.txt 
head prep_for_ref_allele_input.txt 

# chr	start
chr10	100003915

wc -l prep_for_ref_allele_input.txt
# 903,098 prep_for_ref_allele_input.txt (header)
```

Now use R to edit the file further: removing chromosomes 0, XY, MT
```
module load R/3.6.0 
R
library(Rsamtools)
library(BSgenome)
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome")

dat <- read.table("prep_for_ref_allele_input.txt", header=T)

# if you have a genome fasta file, this file can be imported using "FaFile", the format of chromosomes should match that in the fasta file. The labeling has to be the same in the fasta and dat file, e.g. 1, 2, 3 or chr1, chr2, chr3 in both

fasta_file <- FaFile(file='hg19_genome.fa')
# Files with the .fa extension use the FASTA format, which is a file type used for the storage of numerous sequences in one file.

gr1 <- GRanges(dat$chr,IRanges(start=as.numeric(dat$start), end=as.numeric(dat$start)))
# GRanges = genomic ranges that each have a single start and end location on the genome

refbase <- getSeq(fasta_file, gr1)
# getSeq extracts a set of sequences (or subsequences) from a BSgenome or XStringSet object

refbase <- as.data.frame(refbase)$x

dat$REF <- refbase

write.table(dat,file="coding_alleles_ref_allele.txt",quote=F,row.names=F,sep="\t")
q()
n
```

Now merge the new file with the other data
```
# combine column 1 (chr) and 2 (start), separating with : (chr:start).
cut -f 1,2 coding_alleles_ref_allele.txt | sed 's/\t/:/g' > coding_alleles_ref_allelev2.txt

# add other columns to it
paste coding_alleles_ref_allelev2.txt coding_alleles_ref_allele.txt > coding_alleles_ref_allelev3.txt

# writing file with reference alleles that need to be corrected
paste coding_alleles_ref_allelev3.txt exonic_splicing_only_no_dup.txt > to_resolve_REF.txt

wc -l to_resolve_REF.txt
# 903,098 to_resolve_REF.txt
```

Compare if reference allele is in column V4 or V5
```
module load R/3.6.0 
R

library(dplyr)
library(data.table)

to_resolve = fread("to_resolve_REF.txt", sep = "\t", header = T)

to_resolve %>% group_by(V4) %>% tally()
# A tibble: 3 × 2
  V4         n
  <chr>  <int>
1 A     429830
2 C     413349
3 G      59918

to_resolve %>% group_by(V5) %>% tally()
# A tibble: 3 × 2
  V5         n
  <chr>  <int>
1 C      60185
2 G     414589
3 T     428323

to_resolve2 = to_resolve %>% mutate(Comparing_ref = if_else(REF == V4 & REF == V5, "TRUETRUE",
                                                        if_else(REF == V4 & REF != V5, "TRUEFALSE",
                                                            if_else(REF != V4 & REF == V5, "FALSETRUE", 
                                                                if_else(REF != V4 & REF != V5, "FALSEFALSE", "EMPTY")))))

to_resolve2 %>% group_by(Comparing_ref) %>% tally()

# A tibble: 3 × 2
  Comparing_ref      n
  <chr>          <int>
1 FALSEFALSE         9
2 FALSETRUE     452215
3 TRUEFALSE     450873
```

In 9 cases, the reference allele is not available in any of the columns - remove those!

```
to_resolve2 = to_resolve2 %>% filter(Comparing_ref != "FALSEFALSE")

# A tibble: 2 x 2
  Comparing_ref      n
  <chr>          <int>
1 FALSETRUE     452215
2 TRUEFALSE     450873
```

And replace V5 with V4, V4 with REF - for the correct allele order
```
to_resolve3 = to_resolve2
to_resolve3 %>% group_by(Comparing_ref) %>% tally()
# A tibble: 2 x 2
  Comparing_ref      n
  <chr>          <int>
1 FALSETRUE     452215
2 TRUEFALSE     450873


to_resolve3 = to_resolve3 %>% mutate(V5 = ifelse(Comparing_ref == "FALSETRUE", V4, V5))
to_resolve3 = to_resolve3 %>% mutate(V4 = ifelse(Comparing_ref == "FALSETRUE", REF, V4))

to_resolve3 = to_resolve3 %>% mutate(Comparing_ref = if_else(REF == V4 & REF == V5, "TRUETRUE",
                                                            if_else(REF == V4 & REF != V5, "TRUEFALSE",
                                                                    if_else(REF != V4 & REF == V5, "FALSETRUE", 
                                                                            if_else(REF != V4 & REF != V5, "FALSEFALSE", "EMPTY")))))


to_resolve3 %>% group_by(Comparing_ref) %>% tally()
# A tibble: 1 x 2
  Comparing_ref      n
  <chr>          <int>
1 TRUEFALSE     903088

write.table(to_resolve3, "to_resolve_REF_v2.txt", row.names = F, quote = F, sep="\t")
q()
n
```

```
wc -l to_resolve_REF_v2.txt 
# 903,089 to_resolve_REF_v2.txt (header)
```

### Re-annotate

```
# change file into normal layout: 
# Chr, start, start, A1, A2

cut -f 5-9 to_resolve_REF_v2.txt | grep -v "V5" > to_annotate_v2.txt
```

New working file "to_annotate_v2.txt"

```
module load annovar

table_annovar.pl to_annotate_v2.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150,clinvar_20200316 \ 
-operation g,f,f -outfile second_pass_23andme_clinvar -nastring .

# new file second_pass_23andme_clinvar.hg19_multianno.txt
```

Sort out space issues
```
sed -i -e 's/ /_/g' second_pass_23andme_clinvar.hg19_multianno.txt

wc -l second_pass_23andme_clinvar.hg19_multianno.txt
# 903089 second_pass_23andme_clinvar.hg19_multianno.txt (header)
```

## 3. Extracting gene-specific information
Start off with just looking at GBA, Parkinson's and Lewy. Grep the terms and write new files:

Write all info on GBA in separate file
```
grep -w GBA second_pass_23andme_clinvar.hg19_multianno.txt > GBA_only.txt
wc -l GBA_only.txt
# 135 GBA_only.txt (no header)
```

Write all info on Parkinson's into separate file
```
grep -e CLNDISDB -e "arkinson" second_pass_23andme_clinvar.hg19_multianno.txt | \
grep -v "Wolff" > PD_only.txt
wc -l PD_only.txt
# 356 PD_only.txt (header)
```

Write all info on Lewy body dementia
```
grep -e CLNDISDB -e "Lewy" second_pass_23andme_clinvar.hg19_multianno.txt | \
grep -v "Wolff" > Lewy_only.txt
wc -l Lewy_only.txt
# 10 Lewy_only.txt (header)
```

Prioritize variants by removing all that say "Benign" or "Benign/Likely_benign" or "Likely_benign" from "CLNSIG" for the files:
PD_only.txt, GBA_only.txt, and Lewy_only.txt

```
library(dplyr)
library(data.table)

PD_only = fread("PD_only.txt")
dim(PD_only)
# 355 16 

PD_only %>% group_by(CLNSIG) %>% tally()

# CLNSIG                                           n
  <chr>                                        <int>
1 Benign                                          25
2 Benign/Likely_benign                            27
3 Conflicting_interpretations_of_pathogenicity    30
4 Likely_benign                                   19
5 Likely_pathogenic                                1
6 Pathogenic                                      25
7 risk_factor                                      3
8 Uncertain_significance                          87

PD_only = PD_only %>% filter(CLNSIG != "Benign"& CLNSIG != "Benign/Likely_benign" & CLNSIG != "Likely_benign")

dim(PD_only)
# 236  16

PD_only %>% group_by(CLNSIG) %>% tally()

# A tibble: 6 x 2
  CLNSIG                                           n
  <chr>                                        <int>
1 Conflicting_interpretations_of_pathogenicity    49
2 Likely_pathogenic                                2
3 Pathogenic                                      45
4 Pathogenic/Likely_pathogenic                     4
5 risk_factor                                      8
6 Uncertain_significance                         128

PD_only %>% group_by(Gene.refGene) %>% tally() %>% arrange(desc(n))
# A tibble: 30 x 2
   Gene.refGene     n
   <chr>        <int>
 1 LRRK2           75
 2 ATP13A2         26
 3 PINK1           22
 4 PRKN            16
 5 SYNJ1           12
 6 PLA2G6           9
 7 SLC6A3           8
 8 SNCAIP           8
 9 FBXO7            7
10 GBA              7
# … with 20 more rows


write.table(PD_only, file="PD_only.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

```
wc -l PD_only.txt
# 237 PD_only.txt
```

Repeat for GBA_only.txt (n=135)

```
GBA_only = fread("GBA_only.txt")
# this file does not contain headers, so add them using rename()

GBA_only = GBA_only %>% rename("Chr" = V1, "Start" = V2, "End" = V3, "Ref" = V4, "Alt" = V5, "Func.refGene" = V6, "Gene.refGene" = V7, "GeneDetail.refGene" = V8, "ExonicFunc.refGene" = V9, "AAChange.refGene" =V10, "avsnp150" = V11, "CLNALLELE" = V12, "CLNDN" = V13, "CLNDISDB" = V14, "CLNREVSTAT" = V15, "CLNSIG" = V16)

colnames(GBA_only)
 [1] "Chr"                "Start"              "End"               
 [4] "Ref"                "Alt"                "Func.refGene"      
 [7] "Gene.refGene"       "GeneDetail.refGene" "ExonicFunc.refGene"
[10] "AAChange.refGene"   "avsnp150"           "CLNALLELEID"       
[13] "CLNDN"              "CLNDISDB"           "CLNREVSTAT"        
[16] "CLNSIG" 

GBA_only %>% group_by(CLNSIG) %>% tally()

# A tibble: 10 x 2
   CLNSIG                                           n
   <fct>                                        <int>
 1 .                                               82
 2 Benign                                           1
 3 Benign/Likely_benign                             1
 4 Conflicting_interpretations_of_pathogenicity     5
 5 Likely_benign                                    2
 6 Likely_pathogenic                                8
 7 Pathogenic                                      25
 8 Pathogenic/Likely_pathogenic                     4
 9 risk_factor                                      1
10 Uncertain_significance                           6

GBA_only = GBA_only %>% filter(CLNSIG != "Benign"& CLNSIG != "Benign/Likely_benign" & CLNSIG != "Likely_benign")

GBA_only %>% group_by(CLNSIG) %>% tally()
# A tibble: 7 x 2
  CLNSIG                                           n
  <fct>                                        <int>
1 .                                               82
2 Conflicting_interpretations_of_pathogenicity     5
3 Likely_pathogenic                                8
4 Pathogenic                                      25
5 Pathogenic/Likely_pathogenic                     4
6 risk_factor                                      1
7 Uncertain_significance                           6

GBA_only %>% group_by(Gene.refGene) %>% tally()
# A tibble: 1 x 2
  Gene.refGene     n
  <fct>        <int>
1 GBA            131

write.table(GBA_only, file="GBA_only.txt", quote=FALSE,row.names=F,sep="\t")
q()

wc -l GBA_only.txt 
# 132 GBA_only.txt
```

And Lewy_only.txt

```
R
library(data.table)
library(dplyr)

Lewy = read.table("Lewy_only.txt", sep = "\t", header = T)

Lewy %>% group_by(CLNSIG) %>% tally()
# A tibble: 1 x 2
  CLNSIG         n
  <fct>      <int>
1 Pathogenic     6

Lewy %>% group_by(Gene.refGene) %>% tally()
# A tibble: 3 x 2
  Gene.refGene     n
  <fct>        <int>
1 GBA              2
2 SNCA             2
3 SNCB             2

# no filtering needed
wc -l Lewy_only.txt
# 10 Lewy_only.txt
```

### Merge all files, using full join
```
# continue in R

PD_only = fread("PD_only.txt")
dim(PD_only)
# [1] 236  16

GBA_only = fread("GBA_only.txt")
dim(GBA_only)
# [1] 131  16

Lewy_only = fread("Lewy_only.txt")
dim(Lewy_only)
# [1]  9 16

join = full_join(PD_only, GBA_only)
# note that GBA adds another column "CLNALLELE"
dim(join)
# [1] 360  17


join = full_join(join, Lewy_only)
dim(join)
# [1] 363  17
```

```
#identify duplicates that might be introduced because of this extra columns
join %>% group_by(Start) %>% tally() %>% filter(n>1)
# A tibble: 0 x 2
# … with 2 variables: Start <int>, n <int>

join %>% group_by(avsnp150) %>% tally() %>% filter(n>1)
# A tibble: 1 x 2
  avsnp150     n
  <chr>    <int>
1 .           42

## they're all GBA!
```

Write new file
```
write.table(join, file="variants_of_interest.txt", quote=FALSE,row.names=F,sep="\t")
q()
```

```
wc -l variants_of_interest.txt
# 364 variants_of_interest.txt (header) - 363 variants remaining
```

Now extract 23andme annotation format
```
# write columns 5 (scaffold) and 6 (position) into new file and separate them by ":" to match 23andme format 
cut -f 5,6 all_snp_info.txt | sed 's/\t/:/g' > var_name.txt

# write assay.name into new file
cut -f 4 all_snp_info.txt > assay_name.txt

# and combine both files to new file "to_import_R.txt"
paste assay_name.txt var_name.txt > to_import_R.txt
```

Now edit files in R by adding one more column and compare both files
```
module load R/3.6.0 
R
require(data.table)
library(dplyr)
data <- fread("to_import_R.txt",header=T)
voi <- fread("variants_of_interest.txt",header=T)

# add header to second column
voi2 = tidyr::unite(voi, "scaffold:position", c(Chr:Start), sep = ":")

double_check <- merge(data, voi2, by="scaffold:position")
dim(double_check)
# [1] 363  17

write.table(double_check,file="variants_of_interest_with_23andme_ID.txt",quote=F,row.names=F,sep="\t")
q()
n

# file contains 364 variants with 18 columns
```

```
wc -l variants_of_interest_with_23andme_ID.txt
364 variants_of_interest_with_23andme_ID.txt (header)
```

Extract gene-specific info, based on review Blauwendraat et al 2019 (PMID: 31521533)
```
grep -w ATP13A2 second_pass_23andme_clinvar.hg19_multianno.txt > ATP13A2_only.txt
grep -w DNAJC13 second_pass_23andme_clinvar.hg19_multianno.txt > DNAJC13_only.txt
grep -w DNAJC6 second_pass_23andme_clinvar.hg19_multianno.txt > DNAJC6_only.txt
grep -w EIF4G1 second_pass_23andme_clinvar.hg19_multianno.txt > EIF4G1_only.txt
grep -w FBXO7 second_pass_23andme_clinvar.hg19_multianno.txt > FBXO7_only.txt
grep -w GBA second_pass_23andme_clinvar.hg19_multianno.txt > GBA_only.txt
grep -w GIGYF2 second_pass_23andme_clinvar.hg19_multianno.txt > GIGYF2_only.txt
grep -w HTRA2 second_pass_23andme_clinvar.hg19_multianno.txt > HTRA2_only.txt
grep -w LRP10 second_pass_23andme_clinvar.hg19_multianno.txt > LRP10_only.txt
grep -w LRRK2 second_pass_23andme_clinvar.hg19_multianno.txt > LRRK2_only.txt
grep -w PARK7 second_pass_23andme_clinvar.hg19_multianno.txt > PARK7_only.txt
grep -w PINK1 second_pass_23andme_clinvar.hg19_multianno.txt > PINK1_only.txt
grep -w PLA2G6 second_pass_23andme_clinvar.hg19_multianno.txt > PLA2G6_only.txt
grep -w POLG second_pass_23andme_clinvar.hg19_multianno.txt > POLG_only.txt
grep -w PRKN second_pass_23andme_clinvar.hg19_multianno.txt > PRKN_only.txt
grep -w SNCA second_pass_23andme_clinvar.hg19_multianno.txt > SNCA_only.txt
grep -w SYNJ1 second_pass_23andme_clinvar.hg19_multianno.txt > SYNJ1_only.txt
grep -w TMEM230 second_pass_23andme_clinvar.hg19_multianno.txt > TMEM230_only.txt
grep -w UCHL1 second_pass_23andme_clinvar.hg19_multianno.txt > UCHL1_only.txt
grep -w VPS13C second_pass_23andme_clinvar.hg19_multianno.txt > VPS13C_only.txt
grep -w VPS35 second_pass_23andme_clinvar.hg19_multianno.txt > VPS35_only.txt
```

Merge all files
```
# write header only
head -1 second_pass_23andme_clinvar.hg19_multianno.txt > header.txt

# concatenate all files
cat header.txt ATP13A2_only.txt DNAJC13_only.txt DNAJC6_only.txt EIF4G1_only.txt FBXO7_only.txt GBA_only.txt GIGYF2_only.txt HTRA2_only.txt LRP10_only.txt LRRK2_only.txt PARK7_only.txt PINK1_only.txt PLA2G6_only.txt POLG_only.txt PRKN_only.txt SNCA_only.txt SYNJ1_only.txt TMEM230_only.txt UCHL1_only.txt VPS13C_only.txt VPS35_only.txt > PD_genes.txt

wc -l PD_genes.txt
# 1544 PD_genes.txt (header)
```

Annotate again, including allele frequency 
```
# writing columns 1-5 (Chr	Start	End	Ref	Alt) to a new file
cut -f 1-5 PD_genes.txt > to_annotate_PD_genes.txt

# remove header in R
annotate = read.table("to_annotate_PD_genes.txt", header = T, sep = "\t")
names(annotate) <- NULL
write.table(annotate, "to_annotate_PD_genes.txt", quote = F, row.names = F, sep = "\t")

#annotate
table_annovar.pl to_annotate_PD_genes.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150,clinvar_20200316,gnomad211_genome \
-operation g,f,f,f -outfile third_pass_23andme_clinvar_to_annotate_PD_genes -nastring .

# new working file: "third_pass_23andme_clinvar_to_annotate_PD_genes.hg19_multianno.txt"
wc -l third_pass_23andme_clinvar_to_annotate_PD_genes.hg19_multianno.txt 
# 1544 third_pass_23andme_clinvar_to_annotate_PD_genes.hg19_multianno.txt (header)
```

### Remove synonymous variants and "-AS" genes
```
module load R/3.6.0 
R
library(dplyr)
library(tidyr)
anno = read.table("third_pass_23andme_clinvar_to_annotate_PD_genes.hg19_multianno.txt", sep = "\t", header = T, quote="", fill=FALSE)
dim(anno)
# 1543 33

anno %>% group_by(ExonicFunc.refGene) %>% tally()
# A tibble: 5 x 2
  ExonicFunc.refGene     n
  <fct>              <int>
1 .                    144
2 nonsynonymous SNV    947
3 startloss              3
4 stopgain              43
5 synonymous SNV       406

anno = anno %>% filter(ExonicFunc.refGene != "synonymous SNV")
dim(anno)
# 1137 33

anno %>% group_by(ExonicFunc.refGene) %>% tally()
# A tibble: 4 x 2
  ExonicFunc.refGene     n
  <fct>              <int>
1 .                    144
2 nonsynonymous SNV    947
3 startloss              3
4 stopgain              43

# remove genes "-AS" in $Gene.refGene
anno = tidyr::unite(anno, "combo", c(Chr:Start), sep = ":")
anno = anno[!grepl("-AS", anno$Gene.refGene), ]

write.table(anno, "variants_of_interest_additional.txt", quote=F,row.names=F,sep="\t")
# adding the VP because "cannot open file 'variants_of_interest_additional.txt': Permission denied"
q()
n
```

```
wc -l variants_of_interest_additional.txt
# 1015 variants_of_interest_additional.txt
```

Continue editing 23andme annotation and extracting 23andme format
```
# writing columns 5&6 to new file, joining them together by ":" (scaffold:position)
cut -f 5,6 all_snp_info.txt | sed 's/\t/:/g' > var_name.txt

# write column 4 (assay.name) into new file
cut -f 4 all_snp_info.txt > assay_name.txt

#paste them together and write new file
paste assay_name.txt var_name.txt > to_import_R.txt
```

```
module load R
R
require(data.table)
library(dplyr)

data <- fread("to_import_R.txt",header=T)
voi <- fread("variants_of_interest_additional.txt",header=T)

names(data)[2] <- "combo"

# merge here
double_check <- merge(data, voi, by="combo")
dim(double_check)
# 1014 33

write.table(double_check,file="variants_of_interest_with_23andme_ID_second_pass.txt",quote=F,row.names=F,sep="\t")
```

```
wc -l variants_of_interest_with_23andme_ID_second_pass.txt
# 1015 variants_of_interest_with_23andme_ID_second_pass.txt (header)
```

## 4. Final list of variants
```
cat variants_of_interest_with_23andme_ID.txt variants_of_interest_with_23andme_ID_second_pass.txt \
| sort -u | grep -v "assay.name" > FINAL_variants_of_interest_with_23andme_ID.txt
```

```
wc -l FINAL_variants_of_interest_with_23andme_ID.txt 
# 1377 FINAL_variants_of_interest_with_23andme_ID.txt
```

Double-check in R
```
R
library(dplyr)
library(data.table)

### first file
first = fread("variants_of_interest_with_23andme_ID.txt", header = T)
dim(first)
# 363 17
first$CLNALLELEID = as.character(first$CLNALLELEID)

# check variants in here
first %>% group_by(ExonicFunc.refGene) %>% tally()
# A tibble: 4 × 2
  ExonicFunc.refGene     n
  <chr>              <int>
1 .                     11
2 nonsynonymous_SNV    293
3 stopgain              17
4 synonymous_SNV        42

first = first %>% rename("combo" = 'scaffold:position')
first_ID = first %>% select(combo, assay.name, ExonicFunc.refGene)

### second file
second = fread("variants_of_interest_with_23andme_ID_second_pass.txt", header = T)
dim(second)
# 1014 33

# check variants in here
second %>% group_by(ExonicFunc.refGene) %>% tally()
# A tibble: 4 × 2
  ExonicFunc.refGene     n
  <chr>              <int>
1 .                     21
2 nonsynonymous_SNV    947
3 startloss              3
4 stopgain              43

second$ExonicFunc.refGene[second$ExonicFunc.refGene == "nonsynonymous SNV"] <- "nonsynonymous_SNV"

second_ID = second %>% select(combo, assay.name, ExonicFunc.refGene)

IDjoin = full_join(first_ID, second_ID)
dim(IDjoin)
# 1091 3

## now join other columns to this with left_join
final = left_join(IDjoin, second)
final = left_join(final, first)

dim(final)
# 1091 34

## check duplicates
final %>% group_by(combo) %>% tally() %>% filter(n>1)
# none
final %>% group_by(assay.name) %>% tally() %>% filter(n>1)
# none

## final checks for synonymous variants
final %>% group_by(ExonicFunc.refGene) %>% tally()

# A tibble: 5 × 2
  ExonicFunc.refGene     n
  <chr>              <int>
1 .                     27
2 nonsynonymous_SNV    976
3 startloss              3
4 stopgain              43
5 synonymous_SNV        42

final = final %>% filter(ExonicFunc.refGene != "synonymous_SNV")
dim(final)
# 1049 34

## Clinical significance
final %>% group_by(CLNSIG) %>% tally()

# A tibble: 11 × 2
   CLNSIG                                           n
   <chr>                                        <int>
 1 .                                              521
 2 Benign                                          61
 3 Benign/Likely_benign                            27
 4 Conflicting_interpretations_of_pathogenicity    66
 5 Likely_benign                                   48
 6 Likely_pathogenic                               31
 7 Pathogenic                                      86
 8 Pathogenic/Likely_pathogenic                    18
 9 risk_factor                                      7
10 Uncertain_significance                         149
11 NA                                              35


## will leave all those in, would be interesting to see how impactful they are in our analysis

write.table(final, "FINAL_variants_of_interest_with_23andme_ID.txt", row.names =F, sep = "\t", quote = F)
```

```
wc -l FINAL_variants_of_interest_with_23andme_ID.txt
# 1050 FINAL_variants_of_interest_with_23andme_ID.txt
```

DONE
