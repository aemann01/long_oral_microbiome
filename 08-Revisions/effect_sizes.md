Load libraries, data, install effectsizes

```R
library(phyloseq)
load("../02-diversity/master_phyloseq.RData")
library(effectsize)
options(es.use_symbols = TRUE)
```

Effect sizes for taxonomic turnover analyses

```R
library(compositions)
ps.dat
# pull metadata table
map <- sample_data(ps.dat)
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")

map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
counts <- data.frame(table(map.v1v3$studyID_FDI))
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
dim(map.v1v3)

dates <- as.data.frame(matrix(map.v1v3$visit_date, ncol=2, byrow=TRUE))
# average difference between date of visit one and visit two
mean(as.Date(as.character(dates$V2), format="%m/%d/%Y")-as.Date(as.character(dates$V1), format="%m/%d/%Y"))

counts <- otu_table(ps.dat)[row.names(map.v1v3),]
dim(counts)
# remove asvs with zero count after filtering
counts <- counts[,colSums(counts) > 0]
dim(counts)
# make sure data is all numeric
counts <- apply(counts,c(1,2),function(x) as.numeric(as.character(x)))
# perform clr transformation to the data
counts.clr <- clr(counts)
# distance
counts.dist <- dist(counts.clr, method="euclidean")

# compute taxonomic turnover
library(reshape2)
library(dplyr)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE","FALSE")
# get only those that are duplicated samples, clean up again
matches <- counts.melt[counts.melt$match == "TRUE",]
matches <- as.data.frame(matches[,-2])
matches <- matches[!duplicated(matches),]
# get temporary metadata file to merge in with your distance metrics
temp <- as.data.frame(cbind(map.v1v3$studyID_FDI, map.v1v3$hiv_status))
temp <- temp[!duplicated(temp),]
# final dataframe
dist_diff <- merge(x=matches, y=temp, by.x="Var1", by.y="V1")
# some summary statistics on these distances
group_by(dist_diff, V2) %>%
summarise(
count = n(),
median = median(value, na.rm = TRUE),
SD = sd(value, na.rm = TRUE)
)

# split into two groups for cohen's d
sub <- dist_diff %>% filter(V2 %in% c("HI", "HUU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |        95% CI
# -------------------------
# -0.33     | [-0.66, 0.01]
interpret_cohens_d(-0.33, rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)

# HI vs HEU
sub <- dist_diff %>% filter(V2 %in% c("HI", "HEU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |        95% CI
# -------------------------
# 5.86e-03  | [-0.34, 0.35]

interpret_cohens_d(5.86e-03, rules = "cohen1988")
# [1] "very small"
# (Rules: cohen1988)

# HI vs HEU
sub <- dist_diff %>% filter(V2 %in% c("HEU", "HUU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |        95% CI
# -------------------------
# 5.86e-03  | [-0.34, 0.35]

interpret_cohens_d(-0.35, rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)

### Healthy teeth only
ps.dat.h <- subset_samples(ps.dat, tooth_health=="H")
# pull metadata table
map <- sample_data(ps.dat.h)
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
counts <- data.frame(table(map.v1v3$studyID_FDI))
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
counts <- otu_table(ps.dat)[row.names(map.v1v3),]
counts <- counts[,colSums(counts) > 0]
counts <- apply(counts,c(1,2),function(x) as.numeric(as.character(x)))
# perform clr transformation to the data
counts.clr <- clr(counts)
# distance
counts.dist <- dist(counts.clr, method="euclidean")
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE","FALSE")
# get only those that are duplicated samples, clean up again
matches <- counts.melt[counts.melt$match == "TRUE",]
matches <- as.data.frame(matches[,-2])
matches <- matches[!duplicated(matches),]
# get temporary metadata file to merge in with your distance metrics
temp <- as.data.frame(cbind(map.v1v3$studyID_FDI, map.v1v3$hiv_status))
temp <- temp[!duplicated(temp),]
# final dataframe
dist_diff <- merge(x=matches, y=temp, by.x="Var1", by.y="V1")

# split into two groups for cohen's d
sub <- dist_diff %>% filter(V2 %in% c("HI", "HUU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |        95% CI
# -------------------------
# -0.30     | [-0.68, 0.08]
interpret_cohens_d(-0.30 , rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)

# HI vs HEU
sub <- dist_diff %>% filter(V2 %in% c("HI", "HEU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |        95% CI
# -------------------------
# -0.13     | [-0.51, 0.26]

interpret_cohens_d(-0.13 , rules = "cohen1988")
# [1] "very small"
# (Rules: cohen1988)

# HI vs HEU
sub <- dist_diff %>% filter(V2 %in% c("HEU", "HUU"))
cohens_d(value ~ V2, data = sub)
# Cohen's d |         95% CI
# --------------------------
# -0.45     | [-0.83, -0.06]

interpret_cohens_d(-0.45, rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)
```

Effect sizes for CD4 counts between HIV status groups

```R
load("../02-diversity/master_phyloseq.RData")
ps.dat
ps.dat <- subset_samples(ps.dat, cd4_count != "unknown")
map <- as.data.frame(phyloseq::sample_data(ps.dat))
map$hiv_status <- factor(map$hiv_status, levels=c("HI", "HEU", "HUU"))
map <- data.frame(map)

# first HI vs HUU
sub <- map %>% filter(hiv_status %in% c("HI", "HUU"))
cohens_d(as.numeric(cd4_count) ~ hiv_status, data = sub)
# Cohen's d |         95% CI
# --------------------------
# -0.40     | [-0.51, -0.29]

interpret_cohens_d(-0.40, rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)

sub <- map %>% filter(hiv_status %in% c("HI", "HEU"))
cohens_d(as.numeric(cd4_count) ~ hiv_status, data = sub)
# Cohen's d |         95% CI
# --------------------------
# -0.33     | [-0.44, -0.22]

interpret_cohens_d(-0.33 , rules = "cohen1988")
# [1] "small"
# (Rules: cohen1988)

sub <- map %>% filter(hiv_status %in% c("HUU", "HEU"))
cohens_d(as.numeric(cd4_count) ~ hiv_status, data = sub)
# Cohen's d |        95% CI
# -------------------------
# -0.10     | [-0.21, 0.02]

interpret_cohens_d(-0.10 , rules = "cohen1988")
# [1] "very small"
# (Rules: cohen1988)
```