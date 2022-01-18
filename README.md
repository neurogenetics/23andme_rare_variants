## 23andme rare variants

 - **Project:** Assessment of Parkinson's disease pathogenic genes and risk variants
 - **Author(s):** Vanessa and Cornelis
 - **Date Last Updated:** January 2022

---
### Quick Description: 
Multiple genes and rare variants are associated with Parkinson’s disease and have been shown to be either a strong risk factor or to cause disease. However, the majority of genes are under heavy debate about their actual pathogenicity and this analysis aims to assess what genes and rare variants are actually associated with Parkinson’s disease.

### Motivation/Goals:
1) Understanding the data, where is it located, and providing summary.

## Structure of README:

### [1. Understanding data and subset data](#1-understanding-the-data)
This section goes through:
- understanding the data and where is the data located

### [2. Cleaning the data](#2-cleaning-the-data)
How to clean the data and extract relevant information

## 1. Understanding the data
Data can be located in /data/CARD/projects/23andme_annotation. Note: the data is in genome build hg19.
All SNPs are in the file "all_snp_info.txt".

```
#get node
sinteractive --mem=100g --cpus-per-task=20

cd /data/CARD/projects/23andme_annotation
wc -l all_snp_info_txt

# 64523260 all_snp_info.txt
```

The file contains 64,523,260 rows, equaling the same number of SNPs.

```
head -n 1 all_snp_info.txt 
```
File layout: 

| all.data.id | gt.data.id | im.data.id | assay.name | scaffold | position | alleles | ploidy | cytobandgene.context | is.v1 | is.v2 | is.v3 | is.v4 | is.v5 | h550 | omni | strand |
|-------------|------------|------------|------------|----------|----------|---------|--------|----------------------|-------|-------|-------|-------|-------|------|------|--------|


## 2. Cleaning the data
Reduce file to only selected columns and annotate using ANNOVAR
```
cut -f 4,5,6,7,10 all_snp_info.txt > short_all_snp_info.txt
head -n 1 short_all_snp_info.txt 
```
| assay.name | scaffold | position | alleles | gene.context |
|------------|----------|----------|---------|--------------|

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
```

| scaffold | position | position | alleles | assay.name | gene.context |
|----------|----------|----------|---------|------------|--------------|

## 3. Annotation

```
module load annovar

#gene build hg19

table_annovar.pl input_annovar_first_pass.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150 \
-operation g,f -outfile first_pass_23andme -nastring .
```
This step creates the file "first_pass_23andme.hg19_multianno.txt", which looks like this:

| Chr | Start | End | Ref | Alt | Func.refGene | Gene.refGene | GeneDetail.refGene | ExonicFunc.refGene | AAChange.refGene | avsnp150 |
|-----|-------|-----|-----|-----|--------------|--------------|--------------------|--------------------|------------------|----------|

refGene from ANNOVAR adds these columns: Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, while avsnp150 adds rsnumbers using dbSNP annotations. 


```
wc -l first_pass_23andme.hg19_multianno.txt 
# 64,523,261 first_pass_23andme.hg19_multianno.txt
# still same number as before - check
```

Check how many unique rows 
```
sort first_pass_23andme.hg19_multianno.txt | uniq | wc -l
# failed for me but should be 63,825,064
```

Create new file with only unique rows and variant identifier, combining columns Chr and Start, separating them with :

```
cut -f 1,2 first_pass_23andme.hg19_multianno.txt | sed 's/\t/:/g' | sort | uniq -d > multiallelics.txt
wc -l multiallelics.txt
# 671848 multiallelics.txt
```
Some alleles are labelled as D/I but are they useful?
Write them into a new file by replacing D with A and I with G and then annotating again with ANNOVAR:
```
awk '$4 == "D"' input_annovar_first_pass.txt > input_annovar_deletions.txt
wc -l input_annovar_deletions.txt 
# 4,124,008 input_annovar_deletions.txt

#replace the alleles
sed -i -e 's/D/A/g' input_annovar_deletions.txt
sed -i -e 's/I/G/g' input_annovar_deletions.txt
```

```
table_annovar.pl input_annovar_deletions.txt $ANNOVAR_DATA/hg19/ \
-buildver hg19 -protocol refGene,avsnp150 \
-operation g,f -outfile input_annovar_deletions_annovar_AG_converted -nastring .

#new file written input_annovar_deletions_annovar_AG_converted.hg19_multianno.txt
```
From the annotation, it looks like these variants could be useful but without alleles, the data can not be reliably analysed. So further analysis is based on alleles A-C-T-G and ignores D/I.

Extract exonic variants and write into a new file
```
grep "exonic" first_pass_23andme.hg19_multianno.txt > exonic_only.txt
wc -l exonic_only.txt
# 907876 exonic_only.txt
```
Extract splicing variants and write into a new file
```
grep "splicing" first_pass_23andme.hg19_multianno.txt > splicing_only.txt
wc -l splicing_only.txt 
# 11667 splicing_only.txt
```
Now merge exonic_only.txt and splicing_only.txt
```
cat exonic_only.txt splicing_only.txt | sort -u > exonic_splicing_only.txt
```

And add variant name
```
cut -f 1,2 exonic_splicing_only.txt | sed 's/\t/:/g' > var_name.txt
paste exonic_splicing_only.txt var_name.txt > exonic_splicing_only_v2.txt
```

Fixing space issue
```
sed -i 's/ /_/g' exonic_splicing_only_v2.txt
# this removes the spaces in column 9, e.g. nonsynonymous SNV -> nonsynonymous_SNV
```

Remove multi-allelic variants using antijoin().
```
module load R/3.6.0 
R
require(data.table)
library(dplyr)
data <- fread("exonic_splicing_only_v2.txt",header=F)
dim(data)
# 918876 rows, 12 columns
dup <- fread("multiallelics.txt",header=F)
dim(dup)
# 671848 rows, 1 column
names(dup)[1] <- "V12"

data_no_dup <- anti_join(data, dup, by="V12")
dim(data_no_dup)
# 903097 rows, 12 columns
double_check <- merge(data, dup, by="V12")
dim(double_check)
# 15779 rows, 12 columns
```
Only those variants that were in dup *and* in data were removed from data. Since the multiallelics.txt and exonix_splicing_only_v2.txt files are both based on first_pass_23andme.hg19_multianno.txt, the number of rows left in the file makes sense.

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
```

| chr   | start     |
|-------|-----------|
| chr10 | 100003915 |
| chr10 | 100010909 |
| chr10 | 100010922 |
| chr10 | 100011387 |

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
#combine column 1 (chr) and 2 (start), separating with : (chr:start).
cut -f 1,2 coding_alleles_ref_allele.txt | sed 's/\t/:/g' > coding_alleles_ref_allelev2.txt

# add other columns to it
paste coding_alleles_ref_allelev2.txt coding_alleles_ref_allele.txt > coding_alleles_ref_allelev3.txt

# writing file with reference alleles that need to be corrected
paste coding_alleles_ref_allelev3.txt exonic_splicing_only_no_dup.txt > to_resolve_REF.txt
```

Compare if reference allele is in column V4 or V5
```
module load R/3.6.0 
R

to_resolve = read.table("to_resolve_REF.txt", sep = "\t", header = T)

levels(to_resolve$V4)
levels(to_resolve$V5)

to_resolve$V4 = as.character(to_resolve$V4)
to_resolve %>% group_by(V4) %>% tally()
# A tibble: 3 × 2
  V4         n
  <chr>  <int>
1 A     429830
2 C     413349
3 G      59918

to_resolve$V5 = as.character(to_resolve$V5)
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
```

And replace V5 with V4, V4 with REF - for the correct allele order
```
to_resolve3 = to_resolve2
to_resolve3 %>% group_by(V4) %>% tally()

to_resolve3$REF = as.character(to_resolve3$REF)

to_resolve3 = to_resolve3 %>% mutate(V5 = ifelse(Comparing_ref == "FALSETRUE", V4, V5))
to_resolve3 = to_resolve3 %>% mutate(V4 = ifelse(Comparing_ref == "FALSETRUE", REF, V4))

to_resolve3 = to_resolve3 %>% mutate(Comparing_ref = if_else(REF == V4 & REF == V5, "TRUETRUE",
                                                            if_else(REF == V4 & REF != V5, "TRUEFALSE",
                                                                    if_else(REF != V4 & REF == V5, "FALSETRUE", 
                                                                            if_else(REF != V4 & REF != V5, "FALSEFALSE", "EMPTY")))))

to_resolve3 %>% group_by(Comparing_ref) %>% tally()

# A tibble: 1 × 2
  Comparing_ref      n
  <chr>          <int>
1 TRUEFALSE     903088
```







