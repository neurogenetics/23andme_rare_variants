## check how many of the 1091 variants are in the summary stats provided by 23andme

```
cd /data/CARD/projects/23andme_annotation
# files are also sitting here: https://drive.google.com/drive/u/0/folders/1NLzikBOcMyyrl3tQRV6v_SbPrRNaXpqk
```

Add gene name to list of 1091 variants
```
# go back into previous script and work with files "variants_of_interest_with_23andme_ID_refGene.txt" and "variants_of_interest_with_23andme_ID_second_pass_VP.txt"
# replicate some parts but keep refGene column
module load R/3.6.0 
R
require(data.table)
library(dplyr)
data <- fread("to_import_R.txt",header=T)
voi <- fread("variants_of_interest.txt",header=T)

# add header to second column
voi = tidyr::unite(voi, "scaffold:position", c(Chr:Start), sep = ":")

double_check <- merge(data, voi2, by="scaffold:position")
dim(double_check)
[1] 363  18
write.table(double_check,file="variants_of_interest_with_23andme_ID_refGene.txt",quote=F,row.names=F,sep="\t")
q()
n
# file contains 364 variants with 2 columns
```

Now merge this set of variants with second set of variants to again generate final file but with 1 more column
```
variants = fread("variants_of_interest_with_23andme_ID_refGene.txt", header = T)
names(variants)[1] <- "combo"

second_pass = fread("variants_of_interest_with_23andme_ID_second_pass_VP.txt", header = T)
write.table(variants,file="variants_of_interest_with_23andme_ID_refGene.txt",quote=F,row.names=F,sep="\t")

fulljoin = full_join(variants, variants_of_interest_with_23andme_ID_second_pass_VP.txt)
dim(fulljoin)
# 1091 variants

write.table(fulljoin,file="FINAL_variants_of_interest_with_23andme_VP_refGene.txt",quote=F,row.names=F,sep="\t")

q()
n
```






```
FINAL = read.table("FINAL_variants_of_interest_with_23andme_ID.txt", header = F, sep = "\t")
head(FINAL)
FINAL = FINAL %>% rename("assay.name" = V2)

results_23 = read.table("results_for_collaborators.csv", header = T, sep = ",")
# note that the results file is in GRCh38
head(results_23)

check = merge(FINAL, results_23)
check = check %>% select(V1, assay.name) %>% rename(combo = V1)

# include column to indicate if variant from results_for_collaborators.csv was in FINAL_variants_of_interest_with_23andme_ID.txt
check$included = "YES"
Final_new = read.table("List_of_1091variants_refGene_23andme.txt", header = T, sep = "\t")
Final_new = full_join(Final_new, check)
```

Cross-check some variants

```
Final_new %>% filter(grepl("rs1043424", FINAL$assay.name))
results_23 %>% filter(grepl("rs1043424", results_23$assay.name))

Final_new %>% filter(grepl("rs104893877", FINAL$assay.name))
results_23 %>% filter(grepl("rs104893877", results_23$assay.name))

Final_new %>% filter(grepl("rs1064651", FINAL$assay.name))
results_23 %>% filter(grepl("rs1064651", results_23$assay.name))
```

Only 504 out of 1091 variants included by 23andme in summary stats file - why?

Add allele frequency for cases and controls to the file
```
results_23 = results_23 %>% mutate(dose.b.0_freq = 1-(dose.b.0/2),
                                   dose.b.1_freq = 1-(dose.b.1/2))
                                   
# where frequency is >0.5, remove "-1" so it's easier to interpret
results_23 = results_23 %>% mutate(dose.b.0_freq = if_else(dose.b.0_freq > 0.5, (dose.b.0/2), dose.b.0_freq),
                      dose.b.1_freq = if_else(dose.b.1_freq > 0.5, (dose.b.1/2), dose.b.1_freq)) 
                      
# change column order
results_23 = results_23 %>% select(assay.name:dose.b.0, dose.b.0_freq, AA.1:dose.b.1, dose.b.1_freq:p.batch)

#and write new file
write.table(results_23, "results_for_collaborators_freq.csv", row.names = F, quote = F, sep = "\t")
```

Add gene info to results file, so change file layout for annotation
```
results_23_short = results_23 %>% select(scaffold,position, alleles) %>% rename(Chr = scaffold, Start = position, End = position)
results_23_short$Start = results_23$position
results_23_short = results_23_short %>% select(Chr, Start, End, alleles)
results_23_short = separate(results_23_short, alleles, c("Ref", "Alt"), sep = "/")

# file is now in format chr, Start, End, Ref, Alt
# write new file
write.table(results_23_short, "results_for_collaborators_freq_annotate.txt", quote=F, row.names = F, sep = "\t")
q()
n
```

And annotate

```
module load annovar

table_annovar.pl results_for_collaborators_freq_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol refGene,avsnp150 \
-operation g,f -outfile first_pass_results_for_collaborators_freq -nastring .

# which outputs file "first_pass_results_for_collaborators_freq.hg38_multianno.txt"
```

Merge with results_for_collaborators_freq.csv to add Gene.refGene column
```
R
library(dplyr)
library(data.table)
results_23 = read.table("results_for_collaborators.csv", header = T, sep = ",")
to_add = read.table("first_pass_results_for_collaborators_freq.hg38_multianno.txt", header = T, sep = ",")

#figure out columns, only combo and Gene.refGene needed
#add merge steps

#write new file
write.table(fulljoin, "results_for_collaborators_freq_refGene.txt", quote=F, row.names = F, sep = "\t")
q()
n
```

Add src column to "List_of_1091variants_refGene_23andme.txt"
```
R
library(dplyr)
library(tidyr)

results_23 = read.table("results_for_collaborators_freq_refGene.txt", header = T, sep = "\t")
variants = read.table("List_of_1091variants_refGene_23andme.txt", header = T, sep = "\t")

results_23_short = results_23 %>% select(scaffold, position, assay.name, src) %>% rename(Chr = scaffold, Start = position, End = position)
results_23_short = unite(results_23_short, "combo", c(Chr, End), sep = ":")
results_23_short = results_23_short %>% select(assay.name, src)

variants_src = left_join(variants, results_23_short)

write.table(variants_src, "List_of_1091variants_refGene_23andme_src.txt", quote = F, row.names = F, sep = "\t")

q()
n
```

In "scr" column, for all G (genotyped) variants, add gt.rate from "/data/CARD/PD/summary_stats/23andMe/GBA/gt_snp_stat.txt"
In "src" column, for all I (imputed) variant, add avg.rsqr and min.rsqr from "/data/CARD/PD/summary_stats/23andMe/GBA/gt_snp_stat.txt"
```
R
library(dplyr)
library(data.table)

GTrate = fread("/data/CARD/PD/summary_stats/23andMe/GBA/gt_snp_stat.txt", header = T)
avgmin = fread("/data/CARD/PD/summary_stats/23andMe/GBA/im_snp_stat.txt", header = T)
List = fread("List_of_1091variants_refGene_23andme_src.txt", header = T)

GTrate_short = GTrate %>% select(assay.name, gt.rate)
avgmin_short = avgmin %>% select(assay.name, avg.rsqr, min.rsqr)

List_gt = left_join(List, GTrate_short)

List_gt_avgmin = left_join(List_gt, avgmin_short)
List_gt_avgmin = List_gt_avgmin %>% mutate(included = ifelse(is.na(included), "NO", included))

List_gt_avgmin %>% filter(!is.na(gt.rate)| !is.na(avg.rsqr) | !is.na(min.rsqr)) %>% group_by(included) %>% tally()
  included     n
  <chr>    <int>
1 NO         373
2 YES        504
# at least n=373 variants missing from the list of variants provided by 23andme but there's information on them in other files

# add column to indicate these variants 
List_gt_avgmin = List_gt_avgmin %>% mutate(should_be_available = if_else(included == "NO" & !is.na(gt.rate) | included == "NO" & !is.na(avg.rsqr) | included == "NO" & !is.na(min.rsqr), "to provide", if_else(included == "NO" & is.na(gt.rate) & is.na(avg.rsqr) & is.na(min.rsqr), "missing altogether", "already provided")))

List_gt_avgmin %>% group_by(included, should_be_available) %>% tally()
# A tibble: 3 × 3
# Groups:   included [2]
  included should_be_available     n
  <chr>    <chr>               <int>
1 NO       missing altogether    213
2 NO       to provide            373
3 YES      already provided      504

# missing altogether = this variant is missing from summary stats file
# to provide = we have gt.rate, avg.rsqr, min.rsqr information of this variant, so it clearly exists somewhere in their data and they just didn't provide it to us
# already provided = this variant was provided in the summary stats file "results_for_collaborators"

# set all genotyped ones to NA in "gt.rate_Gedit", if not already
List_gt_avgmin %>% filter(src == "G" & is.na(gt.rate)) %>% tally() 
# n 0

# set all imputed ones to NA in "avg.rsqr_Iedit" & "min.rsqr_Iedit", if not already
List_gt_avgmin %>% filter(src == "I" & is.na(avg.rsqr)) %>% tally() 
# n 0
List_gt_avgmin %>% filter(src == "I" & is.na(min.rsqr)) %>% tally() 
# n 0 

write.table(List_gt_avgmin, "List_of_1091variants_refGene_23andme_src_gtavgmin_to_provide.txt", quote = F, row.names = F, sep = "\t")
q()
n
```


