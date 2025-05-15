```bash
cd /home/allie/long_oral_microbiome/07-Revisions-10.7.24
```

Pull any teeth that went from E to D

```R
library(phyloseq, verbose=F)
library(ggplot2, verbose=F)
library(compositions, verbose=F)
library(reshape2, verbose=F)
library(tidyr, verbose=F)
library(dplyr, verbose=F)
library(ggdist, verbose=F)
library(ggpubr, verbose=F)
# load in master phyloseq R data object
load("../02-diversity/master_phyloseq.RData")
# check phyloseq object
ps.dat
# pull metadata table
map <- data.frame(sample_data(ps.dat))

# dataframe for visit one to visit two 
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
# pull visit one and visit two
map.v1v2 <- map[map$visit_num == 1 | map$visit_num == 2,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v1v2$studyID_FDI))
# pull records that occur more than one time
map.v1v2 <- map.v1v2[map.v1v2$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v1v2 <- map.v1v2[order(map.v1v2$studyID_FDI),]
# how many samples meet these criteria?
dim(map.v1v2)

# same thing visit one and three
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
# pull visit one and visit two
map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v1v3$studyID_FDI))
# pull records that occur more than one time
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
# how many samples meet these criteria?
dim(map.v1v3)

# visit two visit three
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
# pull visit one and visit two
map.v2v3 <- map[map$visit_num == 2 | map$visit_num == 3,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v2v3$studyID_FDI))
# pull records that occur more than one time
map.v2v3 <- map.v2v3[map.v2v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v2v3 <- map.v2v3[order(map.v2v3$studyID_FDI),]
# how many samples meet these criteria?
dim(map.v2v3)
```

Now pull any samples that went from E to D at any time point

```R
# keep row names
map.v1v2 <- map.v1v2 %>% mutate(original_row = rownames(map.v1v2))
ed_map.v1v2 <- map.v1v2 %>% filter(tooth_health %in% c("E", "D")) %>% group_by(study_id) %>% filter(n() > 1) %>% ungroup()
# v1 to v3
map.v1v3 <- map.v1v3 %>% mutate(original_row = rownames(map.v1v3))
ed_map.v1v3 <- map.v1v3 %>% filter(tooth_health %in% c("E", "D")) %>% group_by(study_id) %>% filter(n() > 1) %>% ungroup()
# v2 to v3
map.v2v3 <- map.v2v3 %>% mutate(original_row = rownames(map.v2v3))
ed_map.v2v3 <- map.v2v3 %>% filter(tooth_health %in% c("E", "D")) %>% group_by(study_id) %>% filter(n() > 1) %>% ungroup()
# subset sample list
comb_ids <- unique(c(ed_map.v1v2$original_row, ed_map.v1v3$original_row, ed_map.v2v3$original_row))
```

Subset the original phyloseq object and run some analyses on the community changes

```R
ps.dat.sub <- prune_samples(comb_ids, ps.dat)
# clr transformation
ps.dat.clr <- microbiome::transform(ps.dat.sub, transform="clr", target="OTU")
```

Ordinate and get some stats by HIV status group

HI first

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HI")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_health)
pdf("bdiv_hi.EtoD.pdf")
plot_ordination(split, ordcap, type="samples", color="tooth_health") + 
    theme_minimal() +
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()

# calculate significance
library(ecole)
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_health, method="euclidean")

# Now performing 1 pairwise comparisons. Percent progress:
# 100 ...

# p-adjust method: bonferroni

#    pairs SumOfSqs  F.Model         R2 pval p.adj
# 1 D vs E 6015.181 1.139223 0.01210147 0.08  0.08
```

HUU 

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HUU")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_health)
pdf("bdiv_huu.EtoD.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_health") + 
    theme_minimal() +
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()

# calculate significance
library(ecole)
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_health, method="euclidean")

# Now performing 1 pairwise comparisons. Percent progress:
# 100 ...

# p-adjust method: bonferroni

#    pairs SumOfSqs  F.Model         R2  pval p.adj
# 1 E vs D 6003.794 1.155751 0.02058117 0.078 0.078
```

HEU

```R
split <- subset_samples(ps.dat.clr, hiv_status=="HEU")
# CAP aitchison distance
ordcap <- ordinate(split, "CAP", "euclidean", ~tooth_health)
pdf("bdiv_heu.EtoD.pdf")
plot_ordination(split, ordcap, "samples", color="tooth_health") + 
    theme_minimal() +
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()

# calculate significance
library(ecole)
# significantly different?
permanova_pairwise(otu_table(split), grp=sample_data(split)$tooth_health, method="euclidean")

# Now performing 1 pairwise comparisons. Percent progress:
# 100 ...

# p-adjust method: bonferroni

#    pairs SumOfSqs  F.Model        R2  pval p.adj
# 1 D vs E 6097.581 1.087432 0.0271265 0.246 0.246
```

All together with HIV status overlaid

```R
# CAP aitchison distance
ordcap <- ordinate(ps.dat.sub, "CAP", "euclidean", ~hiv_status+tooth_health)
pdf("bdiv_All.EtoD.pdf")
plot_ordination(ps.dat.sub, ordcap, "samples", color="tooth_health", shape="hiv_status") + 
    theme_minimal() +
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
```

Top species in E and D by HIV status group

```R
temp <- tax_glom(ps.dat.sub, taxrank="V8")
ps.rel <- transform_sample_counts(temp, function(x) x/sum(x)*100)

topasvs <- names(sort(taxa_sums(ps.rel), TRUE)[1:5])
ps.50 <- prune_taxa(topasvs, ps.rel)
p <- plot_bar(ps.50, "V8", fill="V8", facet_grid=hiv_status~tooth_health)

pdf("etod_toptax_plot.pdf")
p + geom_bar(aes(color=V8, fill=V8), stat="identity", position="stack") +
	theme_minimal() + 
	theme(axis.text.x=element_blank())
dev.off()

# by sample
ps.melt <- psmelt(ps.50)

# ok so the previous figure with high S mutans in HI E is driven by a single individual -- make boxplots of top five taxa in both 

pdf("etod_boxplot.pdf")
ggplot(ps.melt, aes(x = hiv_status, y = Abundance, fill = V8)) +
  geom_boxplot() +
  facet_wrap(~ tooth_health) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
dev.off()


What taxa dominate E and D in each group?

```R
### HI only
split <- subset_samples(ps.dat.sub, hiv_status=="HI")
x <- as.matrix(otu_table(split))
y <- as.factor(sample_data(split)$tooth_health)

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_EtoD.HI.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_EtoD.HI.pdf")
bal$`predictions plot`
dev.off()

### HEU only
split <- subset_samples(ps.dat.sub, hiv_status=="HEU")
x <- as.matrix(otu_table(split))
y <- as.factor(sample_data(split)$tooth_health)

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_EtoD.HEU.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_EtoD.HEU.pdf")
bal$`predictions plot`
dev.off()

### HUU only
split <- subset_samples(ps.dat.sub, hiv_status=="HUU")
x <- as.matrix(otu_table(split))
y <- as.factor(sample_data(split)$tooth_health)

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_EtoD.HUU.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_EtoD.HUU.pdf")
bal$`predictions plot`
dev.off()
```

Across all groups?

What bacteria dominate E vs D across all samples?

```R
library(coda4microbiome)
x <- as.matrix(otu_table(ps.dat.sub))
y <- as.factor(sample_data(ps.dat.sub)$tooth_health)

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_EtoD.All.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_EtoD.All.pdf")
bal$`predictions plot`
dev.off()
