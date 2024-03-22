# Spatial stability of the oral microbiome in children with HIV

Previous work by coauthor CA on a subset of samples found that there is a graident of taxa from the back to the front of the mouth with incisors and molars having very distinct communities in H-CF teeth in HUU and HEU but not in HI. These analyses expand on these preliminary results.

## 1. Install required libraries

```R
# install.packages("remotes")
# remotes::install_github("phytomosaic/ecole")
``` 

## 2. Load required libraries

```R
library(ecole)
library(phyloseq)
library(microbiome)
library(tidyr)
```

## 3. Load in phyloseq object created for diversity analyses

```R
# make folder for figures
system("mkdir img")
load("../02-diversity/master_phyloseq.RData")
```

## 4. Filter phyloseq object 

Want to only look at: healthy teeth from individuals with no caries experience, only adult teeth

```R
set.seed(48793)
# replace NAs with unknown
sample_data(ps.dat)$tooth_type <- sample_data(ps.dat)$tooth_type %>% replace_na("unknown")
# remove samples without tooth type info
ps.dat <- subset_samples(ps.dat, tooth_type!="unknown")
# only looking at adult teeth
ps.dat <- subset_samples(ps.dat, tooth_age=="adult")
# rename premolars as molars to help with low sample size
sample_data(ps.dat)$tooth_type[sample_data(ps.dat)$tooth_type == "premolar"] <- "molar"
ps.dat
```

## 5. CLR transformation

```R
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
```

## 6. Incisors vs molars in H-CF teeth only

```R
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
```

HI samples

```R
split <- subset_samples(temp, hiv_status=="HI")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```

HEU samples

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HEU")
split <- subset_samples(split, aliquot_type=="H-CF")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```

HUU samples

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HUU")
split <- subset_samples(split, aliquot_type=="H-CF")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```

## 6. Incisors vs molars in H-CF teeth only

```R
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
```

HI samples

```R
split <- subset_samples(temp, hiv_status=="HI")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```

HEU samples

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HEU")
split <- subset_samples(split, aliquot_type=="H-CF")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```

HUU samples

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HUU")
split <- subset_samples(split, aliquot_type=="H-CF")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_type)
pdf("img/bdiv_hi.hcf.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
dev.off()
plot_ordination(split, ordcap, "samples", color="tooth_type") + theme_minimal()
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_type, method="euclidean")
```