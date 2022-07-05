# Generate list of frequencies only for 23andMe file for power calculation

```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
```

```
module load R
R
library(tidyverse)


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

final %>% filter(is.na(pvalue)) %>% tally()
# 151

final = final %>% filter(!is.na(pvalue)) 

write.table(final, "23andMe_variants_only_significant.txt", row.names = F, quote = F, sep = "\t")

## edit file for annotation
data = read.table("23andMe_variants_only_significant.txt", header = T, sep = "\t")
data = data %>% select(CHR.BP.REF.ALT) %>% separate(CHR.BP.REF.ALT, c("CHR", "BP", "REF", "ALT"), sep = ":")
data$End = data$BP
data = data %>% select(CHR, BP, End, REF, ALT)
dim(data)
# 812 5

names(data) <- NULL
write.table(data, "23andMe_812_significant_to_annotate.txt", row.names = F, sep = "\t", quote = F)
```

### Then annotate using gnomad_genome for the frequencies (MAF)
```
module load annovar

table_annovar.pl 23andMe_812_significant_to_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol gnomad_genome \
-operation f -outfile 23andMe_812variants_MAF -nastring .
# new file: 23andMe_812variants_MAF.hg38_multianno.txt
```

### Edit frequency in 23andMe file - it's messy!
```
R
library(tidyverse)
andme = read.table("23andMe_variants_only_significant.txt", header = T, sep = "\t")

# overall how many are genotyped vs imputed?
andme %>% group_by(src) %>% tally()
# A tibble: 2 × 2
  src       n
  <chr> <int>
1 G       303
2 I       509

# how many missing freq.b?
andme %>% filter(is.na(freq.b)) %>% tally()
    n
1 374

# correct those imputed ones with missing freq.b
andme$MAF_23 = andme$freq.b
andme_edit = andme %>% mutate(MAF_23 = if_else(src == "I" & is.na(freq.b), dose.b, MAF_23))
andme_edit = andme_edit %>% select(CHR.BP.REF.ALT, src, dose.b, freq.b, MAF_23,effect)

andme_edit %>% filter(is.na(MAF_23)) %>% tally()
   n
1 0

write.table(andme_edit, "23andMe_edited_freq.txt", row.names = F, sep = "\t", quote=F)
```

### Now join both MAF columns to one file
```
R
library(tidyverse)
data = read.table("23andMe_812variants_MAF.hg38_multianno.txt", header = T, sep = "\t")
data = data %>% select(Chr:Alt, gnomAD_genome_NFE)
data = data %>% unite("CHR.BP.REF.ALT", c(Chr, Start, Ref, Alt), sep = ":")
data = data %>% select(CHR.BP.REF.ALT, gnomAD_genome_NFE) %>% rename("MAF_NFE" = gnomAD_genome_NFE)
dim(data)
# 812 2

# how many have missing values
data %>% filter(MAF_NFE == ".") %>% tally()
    n
1 252

andme = read.table("23andMe_edited_freq.txt", header = T, sep = "\t")
dim(andme)
# 812 2

# how many have missing values
andme %>% filter(is.na(MAF_23)) %>% tally()
    n
1 0

merge = merge(data, andme)
dim(merge)
# 812 3

## FLip alleles, i.e. those >0.5, apply 1-freq
merge = merge %>% mutate(MAF_23_flip = if_else(MAF_23 > 0.50, 1-MAF_23, MAF_23))
merge$MAF_NFE = as.character(merge$MAF_NFE)
merge$MAF_NFE = as.numeric(merge$MAF_NFE)
merge = merge %>% mutate(MAF_NFE_flip = if_else(MAF_NFE > 0.50, 1-MAF_NFE, MAF_NFE))

merge = merge %>% select(CHR.BP.REF.ALT, MAF_23_flip, MAF_NFE_flip, effect) %>% rename("MAF_23" = MAF_23_flip, "MAF_NFE" = MAF_NFE_flip, "Effect_23andMe" = effect)
```

### Add meta-analysis effect column
```
meta = read.table("/data/CARD/projects/23andme_annotation/variants_in_internal_datasets/Results/META_AMP_UKB_23andme_clinvar.txt", header = T, sep = "\t")
meta = meta %>% select(CHR.BP.REF.ALT, Effect) %>% rename("Effect_meta" = Effect)

join = left_join(merge, meta)
dim(join)
# 812 5

join = join %>% mutate(Effect_23andMe_abs = abs(Effect_23andMe), Effect_meta_abs = abs(Effect_meta))
write.table(join, "23andMe_812variants_MAF_NFE_effect.txt", row.names = F, sep = "\t", quote=F)
```

### Plot histogram of 23andMe effect (beta)
```
hist(test_data$Effect_23andMe_abs)
```

### Plot frequency from gnomad (NFE) vs 23andMe 
note that despite it saying MAF, these aren't MAFs but only frequencies
```
MAF_23_flip = ggplot(test_data, aes(MAF_NFE, MAF_23)) +
  geom_point() 
ggsave("23andMe_MAF_23_flip_vs_NFE_flip.png", MAF_23_flip, width = 12, height = 5, dpi=300, units = "in")
```

The actual calculations in GAS can be followed here: 
https://hackmd.io/27l0s1BSQZ6e00WewIRNxA?view

In addition to that:

## Derive number of variants with those frequencies:

```
data %>% mutate(Freq_groups = ifelse(MAF_NFE <0.0001, 0.0001,
                                           ifelse(MAF_NFE >= 0.0001 & MAF_NFE <0.0005, 0.0005,
                                                   ifelse(MAF_NFE>=0.0005 & MAF_NFE < 0.001, 0.001,
                                                           ifelse(MAF_NFE>=0.001 & MAF_NFE < 0.005, 0.005,
                                                                   ifelse(MAF_NFE>=0.005 & MAF_NFE <0.01, 0.01,"")))))) %>% group_by(Freq_groups) %>% tally()

```



| OR/MAF             | 0.0001 | 0.0005 | 0.001 | 0.005 | 0.01 |
|--------------------|--------|--------|-------|-------|------|
| Number of variants, n | 252    | 130    | 51    | 63    | 12   |
| 1.5                | 15.5%  | 55.2%  | 84%   | 100%  | 100% |
| 2                  | 38.1%  | 95.9%  | 99.9% | 100%  | 100% |
| 3                  | 79%    | 100%   | 100%  | 100%  | 100% |
| 5                  | 99.2%  | 100%   | 100%  | 100%  | 100% |

## Add actual case numbers
number of cases + 1/4*proxy case = total number of cases

23andMe: 25,034 cases (no proxies) and 3,065,473 controls

0.0001: 25,034*2*0.0001= 5.01 = 5
0.0005: 25,034*2*00005= 25.03 = 25
0.001: 25,034*2*0.001= 50.07 = 50
0.005: 25,034*2*0.005= 250.34 = 250
0.01: 25,034*2*0.01= 500.68 = 501
                                                    

## Final table
This table now includes number of variants with different frequencies, number of cases (from a total of ~25k) among those variants with different frequencies, and power calculations for different OR and MAF scenarios. Based on 23andMe data only.

| OR/MAF                  | 0.0001  | 0.0005 | 0.001  | 0.005  | 0.01  |
|-------------------------|---------|--------|--------|--------|-------|
| Number of variants      | 252     | 130    | 51     | 63     | 12    |
| Number of cases per MAF | 5 | 25 | 50 | 250 | 501 |
| 1.5                     | 15.5%   | 55.2%  | 84%    | 100%   | 100%  |
| 2                       | 38.1%   | 95.9%  | 99.9%  | 100%   | 100%  |
| 3                       | 79%     | 100%   | 100%   | 100%   | 100%  |
| 5                       | 99.2%   | 100%   | 100%   | 100%   | 100%  |

We're reaching 80% power at OR = 3 for all MAFs!


## Calculate with genpwr in R

```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets/Power_calculations/
```

## Power calculation per variant
```
### as suggested by Karl 
# /data/CARD/projects/23andme_annotation/variants_in_internal_datasets/Power_calculations

library(tidyverse)
library(dplyr)
library(ggplot2)
library(devtools)
library(genpwr)

## Taking all cases from meta analysis, accounting for proxies which count 1/4 of a case:
# n_cases = 31156.25 = 31156
# n_controls = 3,108,134
# BUT using the conservative approach as suggested by Mike:
# n_cases = 25,034
# n_controls = 100,000

## Read in data and write a paste loop to get the power estimate for every variant
andme = read.table("23andMe_812variants_MAF_NFE_effect.txt", header = T, sep = "\t")
head(andme)

CHR_MAFs = andme %>% select(MAF_23)
class(CHR_MAFs$MAF_23)
```

Calculate power for OR2
```
CHR_MAFs_OR2 = CHR_MAFs %>% filter(MAF_23 > 1e-9)
CHR_MAFs_OR2 = list(CHR_MAFs_OR2)
# results in 809 MAFs

final = data.frame()

for (i in CHR_MAFs_OR2) {
  power = genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                      N=125000, Case.Rate=0.2, k=NULL,
                      MAF=i[[1]], OR=c(2),Alpha=0.05,
                      True.Model=c("Additive"), 
                      Test.Model=c("Additive"))
  #power1 = c(cbind(MAF = power[1, 3], Power_at_Alpha_0.05 = power[1,9]))
  final_OR2 = rbind(final, power)
}

tobindOR2 = andme %>% filter(MAF_23 > 1e-9) %>% select(CHR.BP.REF.ALT, MAF_23, OR_23andMe)
final2 = cbind(tobindOR2,final_OR2)
head(final2)
```

calculate power for OR3
```
CHR_MAFs_OR3 = CHR_MAFs %>% filter(MAF_23 > 1e-4)
CHR_MAFs_OR3 = list(CHR_MAFs_OR3) 
# results in 335 MAFs

final = data.frame()

for (i in CHR_MAFs_OR3) {
  power = genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                      N=125000, Case.Rate=0.2, k=NULL,
                      MAF=i[[1]], OR=c(3),Alpha=0.05,
                      True.Model=c("Additive"), 
                      Test.Model=c("Additive"))
  #power1 = c(cbind(MAF = power[1, 3], Power_at_Alpha_0.05 = power[1,9]))
  final_OR3 = rbind(final, power)
}

tobindOR3 = andme %>% filter(MAF_23 > 1e-4) %>% select(CHR.BP.REF.ALT, MAF_23, OR_23andMe)
final3 = cbind(tobindOR3,final_OR3)
head(final3)
```
Join them together
```
final = left_join(final2, final3)
final = final %>% select(-c(MAF))
```

Add clinical significance to differentiate
```
Annotations = read.table("Edited_AAChange_Metaanalysis_clean_withannotations.txt", header = T, sep = "\t")
Annotations = Annotations %>% select(MarkerName, CLNSIG) %>% rename("CHR.BP.REF.ALT" = MarkerName)
head(Annotations)

# join annotation on the file
powerFile = left_join(final, Annotations)
head(powerFile)

# change power to %
powerFile = powerFile %>% mutate(Power_at_OR2 = Power_at_OR2*100,
                                 Power_at_OR3 = Power_at_OR3*100)

# write table
write.table(powerFile, "809_variants_with_power_at_OR2_and_OR3.txt", row.names = F, sep = "\t", quote = F)

## how many pathogenic or pathogenic/likely pathogenic ones?
powerFile %>% filter(CLNSIG == "Pathogenic" | CLNSIG == "Pathogenic/Likely_pathogenic") %>% tally()
# 71
```

Write one big file with meta-analysis and power calculations
```
## read in meta analysis output with edited variant names
Meta =  read.table("~/Documents/Projects/23andme_rare_variants/FINAL/Meta_analysis/Edited_AAChange_Metaanalysis_clean.txt", header = T, sep = "\t")
head(Meta)
Meta = Meta %>% rename("CHR.BP.REF.ALT" = MarkerName)

Joined = left_join(Meta, powerFile)
head(Joined)

write.table(Joined, "833_variants_final_metafile_with_powers.txt", row.names = F, sep = "\t", quote = F)
```

*mic drop*
