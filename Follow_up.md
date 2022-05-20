## Generate list of frequencies only for 23andMe file for power calculation
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

```R
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

