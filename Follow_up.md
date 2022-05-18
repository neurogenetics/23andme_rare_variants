## Generate frequency files for power calculation

### generate list of MAF only for 23andMe file for power calculation
```
cd /data/CARD/projects/23andme_annotation/variants_in_internal_datasets
```

```
module load R
R
library(tidyverse)

data = read.table("results_collaborators_all_963.txt", header = T, sep = "\t")

## edit file for annotation
data = data %>% select(CHR.BP.REF.ALT) %>% separate(CHR.BP.REF.ALT, c("CHR", "BP", "REF", "ALT"), sep = ":")
data$End = data$BP
data = data %>% select(CHR, BP, End, REF, ALT)
names(data) <- NULL
write.table(data, "963_original_to_annotate.txt", row.names = F, sep = "\t", quote = F)
```

### then annotate using gnomad_genome for the frequencies (not MAF)
```
module load annovar

table_annovar.pl 963_original_to_annotate.txt $ANNOVAR_DATA/hg38/ \
-buildver hg38 -protocol gnomad_genome \
-operation f -outfile 23andMe_963variants_MAF -nastring .
# new file: 23andMe_963variants_MAF.hg38_multianno.txt
```

### now join both frequency columns to one file
```
R
library(tidyverse)
data = read.table("23andMe_963variants_MAF.hg38_multianno.txt", header = T, sep = "\t")
data = data %>% select(Chr:Alt, gnomAD_genome_NFE)
data = data %>% unite("CHR.BP.REF.ALT", c(Chr, Start, Ref, Alt), sep = ":")
data = data %>% select(CHR.BP.REF.ALT, gnomAD_genome_NFE) %>% rename("MAF_NFE" = gnomAD_genome_NFE)
dim(data)
963 2

andme = read.table("results_collaborators_all_963.txt", header = T, sep = "\t")
dim(andme)
963 2

merge = merge(data, andme)
dim(merge)
963 3
```


### FLip alleles, i.e. those >0.5, apply 1-freq
```
merge = merge %>% mutate(freq.b_flip = if_else(freq.b > 0.50, 1-freq.b, freq.b))
merge$MAF_NFE = as.character(merge$MAF_NFE)
merge$MAF_NFE = as.numeric(merge$MAF_NFE)
merge = merge %>% mutate(MAF_NFE_flip = if_else(MAF_NFE > 0.50, 1-MAF_NFE, MAF_NFE))

merge = merge %>% select(CHR.BP.REF.ALT, freq.b_flip, MAF_NFE_flip) %>% rename("MAF_23" = freq.b_flip, "MAF_NFE" = MAF_NFE_flip)
write.table(merge, "23andMe_963variants_MAF_NFE.txt", row.names = F, sep = "\t", quote=F)
```
### Plot to double check
```
Alleles = ggplot(meta, aes(MinFreq, MaxFreq)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggsave("Allele_issues_meta.png", Alleles, width = 12, height = 5, dpi=300, units = "in")
```

