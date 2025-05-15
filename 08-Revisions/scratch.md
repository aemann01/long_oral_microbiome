1. pull teeth that have multiple visits -- v1 vs v2 v2 vs v3 v1 vs v3
since these are not done on the same day, will have to average the time between visits (we have that information)

```R
# load in phyloseq
library(phyloseq)
library(ggplot2)
# load in master phyloseq data
load("../02-diversity/master_phyloseq.RData")
system("mkdir img")
ps.dat
# pull metadata table
map <- sample_data(ps.dat)

# create a column that combines the sample ID and FDI code
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")

# pull visit one and visit two 
map.v1v2 <- map[map$visit_num == 1 | map$visit_num == 2,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v1v2$studyID_FDI))
# pull records that occur more than one time
map.v1v2 <- map.v1v2[map.v1v2$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v1v2 <- map.v1v2[order(map.v1v2$studyID_FDI),]
# sanity check
dim(map.v1v2)
# head(map.v1v2)

# get average number of days between visits
# first melt visit date column into a data frame of two columns
dates <- as.data.frame(matrix(map.v1v2$visit_date, ncol=2, byrow=TRUE))
# average difference between date of visit one and visit two
mean(as.Date(as.character(dates$V2), format="%m/%d/%Y")-as.Date(as.character(dates$V1), format="%m/%d/%Y"))

## visit two and three? might be able to merge them with v1 and v2 if similar time between sampling
map.v2v3 <- map[map$visit_num == 2 | map$visit_num == 3,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v2v3$studyID_FDI))
# pull records that occur more than one time
map.v2v3 <- map.v2v3[map.v2v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v2v3 <- map.v2v3[order(map.v2v3$studyID_FDI),]
# sanity check
dim(map.v2v3)
# head(map.v2v3)

# get average number of days between visits
# first melt visit date column into a data frame of two columns
dates <- as.data.frame(matrix(map.v2v3$visit_date, ncol=2, byrow=TRUE))
# average difference between date of visit one and visit two
mean(as.Date(as.character(dates$V2), format="%m/%d/%Y")-as.Date(as.character(dates$V1), format="%m/%d/%Y"))
# a little longer than visit one to two, but still half as long as visit one to three. keep separate for now, might merge in future?

## same for visit one and three
map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v1v3$studyID_FDI))
# pull records that occur more than one time
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
# reorder
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
# sanity check
dim(map.v1v3)
# head(map.v1v3)

# get average number of days between visits
# first melt visit date column into a data frame of two columns
dates <- as.data.frame(matrix(map.v1v3$visit_date, ncol=2, byrow=TRUE))
# average difference between date of visit one and visit two
mean(as.Date(as.character(dates$V2), format="%m/%d/%Y")-as.Date(as.character(dates$V1), format="%m/%d/%Y"))
```

now we need to pull these data sets from the OTU table, perform a CLR transformation and get a distance matrix. this is following tutorial by thomazbastiaanssen (Tjazi)

```R
library(compositions)
# first pull the v1 and v2 samples from our phyloseq ASV abundance table and save as new object
counts <- otu_table(ps.dat)[row.names(map.v1v2),]
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
# pca 
counts.pca <- prcomp(counts.dist)
# extract the amount of variance for first four components
pc1 <- round(counts.pca$sdev[1]^2/sum(counts.pca$sdev^2),4) * 100
pc2 <- round(counts.pca$sdev[2]^2/sum(counts.pca$sdev^2),4) * 100
pc3 <- round(counts.pca$sdev[3]^2/sum(counts.pca$sdev^2),4) * 100
pc4 <- round(counts.pca$sdev[4]^2/sum(counts.pca$sdev^2),4) * 100

# extract scores for each sample for first four components
pca <- data.frame(PC1 = counts.pca$x[,1], 
                  PC2 = counts.pca$x[,2], 
                  PC3 = counts.pca$x[,3], 
                  PC4 = counts.pca$x[,4])

# add metadata to the pca scores
pca$studyID_FDI <- map.v1v2$studyID_FDI # individual tooth identifier
pca$hiv_status <- map.v1v2$hiv_status # cohort to split frames
pca$tooth_health <- map.v1v2$tooth_health # treatment
pca$visit_num <- map.v1v2$visit_num # time point data

# now can plot the first two components of the PCA
## NOTE: come back to these and clean up!
pdf("img/test_volplot.pdf")
ggplot(pca, aes(x = PC1, 
                y = PC2, 
                fill = hiv_status,
                colour = hiv_status,
                shape = visit_num, 
                group = studyID_FDI)) +  
	geom_line() +
	geom_point(size = 3, col = "black") + 
	facet_wrap(~tooth_health, scales = "free_x", strip.position = "top") +
	scale_fill_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	scale_colour_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	theme_minimal() +
	xlab(paste("PC1: ", pc1,  "%", sep="")) + 
	ylab(paste("PC2: ", pc2,  "%", sep="")) + 
	theme(text = element_text(size = 12)) +
	guides(fill = guide_legend(override.aes = list(shape = 22)))
dev.off()

# compute volatility -- euclidean distance over CLR transformed count data (aitchison distance between pairs)
# first need to melt our distance matrix to get pairs of samples and their distance
library(reshape2)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE", "FALSE")
# get only those that are duplicated samples, clean up again
matches <- counts.melt[counts.melt$match == "TRUE",]
matches <- as.data.frame(matches[,-2])
matches <- matches[!duplicated(matches),]

# get temporary metadata file to merge in with your distance metrics
temp <- as.data.frame(cbind(map.v1v2$studyID_FDI, map.v1v2$hiv_status))
temp <- temp[!duplicated(temp),]
# final dataframe
dist_diff <- merge(x=matches, y=temp, by.x="Var1", by.y="V1")

# some summary statistics on these distances
library(dplyr)
group_by(dist_diff, V2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
)

# are they different?
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HEU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HEU",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)

# boxplots
pdf("img/test_bplot.pdf")
ggplot(dist_diff, aes(x=V2, y=value)) + geom_boxplot() + geom_jitter()
dev.off()
```

What about visit one to visit three?

```R
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
# pca 
counts.pca <- prcomp(counts.dist)
# extract the amount of variance for first four components
pc1 <- round(counts.pca$sdev[1]^2/sum(counts.pca$sdev^2),4) * 100
pc2 <- round(counts.pca$sdev[2]^2/sum(counts.pca$sdev^2),4) * 100
pc3 <- round(counts.pca$sdev[3]^2/sum(counts.pca$sdev^2),4) * 100
pc4 <- round(counts.pca$sdev[4]^2/sum(counts.pca$sdev^2),4) * 100

# extract scores for each sample for first four components
pca <- data.frame(PC1 = counts.pca$x[,1], 
                  PC2 = counts.pca$x[,2], 
                  PC3 = counts.pca$x[,3], 
                  PC4 = counts.pca$x[,4])

# add metadata to the pca scores
pca$studyID_FDI <- map.v1v3$studyID_FDI # individual tooth identifier
pca$hiv_status <- map.v1v3$hiv_status # cohort to split frames
pca$tooth_health <- map.v1v3$tooth_health # treatment
pca$visit_num <- map.v1v3$visit_num # time point data

# now can plot the first two components of the PCA
## NOTE: come back to these and clean up!
pdf("img/test_volplot.pdf")
ggplot(pca, aes(x = PC1, 
                y = PC2, 
                fill = hiv_status,
                colour = hiv_status,
                shape = visit_num, 
                group = studyID_FDI)) +  
	geom_line() +
	geom_point(size = 3, col = "black") + 
	facet_wrap(~tooth_health, scales = "free_x", strip.position = "top") +
	scale_fill_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	scale_colour_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	theme_minimal() +
	xlab(paste("PC1: ", pc1,  "%", sep="")) + 
	ylab(paste("PC2: ", pc2,  "%", sep="")) + 
	theme(text = element_text(size = 12)) +
	guides(fill = guide_legend(override.aes = list(shape = 22)))
dev.off()

# compute volatility -- euclidean distance over CLR transformed count data (aitchison distance between pairs)
# first need to melt our distance matrix to get pairs of samples and their distance
library(reshape2)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE", "FALSE")
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
library(dplyr)
group_by(dist_diff, V2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
)

# are they different?
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HEU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HEU",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)

# boxplots
pdf("img/test_bplot.pdf")
ggplot(dist_diff, aes(x=V2, y=value)) + geom_boxplot() + geom_jitter()
dev.off()

```

What about visit two to visit three?

```R
counts <- otu_table(ps.dat)[row.names(map.v2v3),]
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
# pca 
counts.pca <- prcomp(counts.dist)
# extract the amount of variance for first four components
pc1 <- round(counts.pca$sdev[1]^2/sum(counts.pca$sdev^2),4) * 100
pc2 <- round(counts.pca$sdev[2]^2/sum(counts.pca$sdev^2),4) * 100
pc3 <- round(counts.pca$sdev[3]^2/sum(counts.pca$sdev^2),4) * 100
pc4 <- round(counts.pca$sdev[4]^2/sum(counts.pca$sdev^2),4) * 100

# extract scores for each sample for first four components
pca <- data.frame(PC1 = counts.pca$x[,1], 
                  PC2 = counts.pca$x[,2], 
                  PC3 = counts.pca$x[,3], 
                  PC4 = counts.pca$x[,4])

# add metadata to the pca scores
pca$studyID_FDI <- map.v1v3$studyID_FDI # individual tooth identifier
pca$hiv_status <- map.v1v3$hiv_status # cohort to split frames
pca$tooth_health <- map.v1v3$tooth_health # treatment
pca$visit_num <- map.v1v3$visit_num # time point data

# now can plot the first two components of the PCA
## NOTE: come back to these and clean up!
pdf("img/test_volplot.pdf")
ggplot(pca, aes(x = PC1, 
                y = PC2, 
                fill = hiv_status,
                colour = hiv_status,
                shape = visit_num, 
                group = studyID_FDI)) +  
	geom_line() +
	geom_point(size = 3, col = "black") + 
	facet_wrap(~tooth_health, scales = "free_x", strip.position = "top") +
	scale_fill_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	scale_colour_manual(values=c("HI"="purple", "HEU"="blue", "HUU"="green")) +
	theme_minimal() +
	xlab(paste("PC1: ", pc1,  "%", sep="")) + 
	ylab(paste("PC2: ", pc2,  "%", sep="")) + 
	theme(text = element_text(size = 12)) +
	guides(fill = guide_legend(override.aes = list(shape = 22)))
dev.off()

# compute volatility -- euclidean distance over CLR transformed count data (aitchison distance between pairs)
# first need to melt our distance matrix to get pairs of samples and their distance
library(reshape2)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE", "FALSE")
# get only those that are duplicated samples, clean up again
matches <- counts.melt[counts.melt$match == "TRUE",]
matches <- as.data.frame(matches[,-2])
matches <- matches[!duplicated(matches),]

# get temporary metadata file to merge in with your distance metrics
temp <- as.data.frame(cbind(map.v2v3$studyID_FDI, map.v2v3$hiv_status))
temp <- temp[!duplicated(temp),]
# final dataframe
dist_diff <- merge(x=matches, y=temp, by.x="Var1", by.y="V1")

# some summary statistics on these distances
library(dplyr)
group_by(dist_diff, V2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
)

# are they different?
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HEU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HEU",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)

# boxplots
pdf("img/test_bplot.pdf")
ggplot(dist_diff, aes(x=V2, y=value)) + geom_boxplot() + geom_jitter()
dev.off()
```

So the previous code calculates volatility across all samples but really we want to look at the community in healthy teeth -- is it more or less volatile in kids with HIV?

```R
# first filter so that we are only looking at H-CF samples in our phyloseq object
ps.dat.hcf <- subset_samples(ps.dat, aliquot_type=="H-CF")
# pull metadata table
map <- sample_data(ps.dat.hcf)
# this should look similar to the processing steps above. Basically want to create a new ASV table with duplicate samples across visit one and visit three
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
# pca 
counts.pca <- prcomp(counts.dist)
# extract the amount of variance for first four components
pc1 <- round(counts.pca$sdev[1]^2/sum(counts.pca$sdev^2),4) * 100
pc2 <- round(counts.pca$sdev[2]^2/sum(counts.pca$sdev^2),4) * 100
pc3 <- round(counts.pca$sdev[3]^2/sum(counts.pca$sdev^2),4) * 100
pc4 <- round(counts.pca$sdev[4]^2/sum(counts.pca$sdev^2),4) * 100

# extract scores for each sample for first four components
pca <- data.frame(PC1 = counts.pca$x[,1], 
                  PC2 = counts.pca$x[,2], 
                  PC3 = counts.pca$x[,3], 
                  PC4 = counts.pca$x[,4])

# add metadata to the pca scores
pca$studyID_FDI <- map.v1v3$studyID_FDI # individual tooth identifier
pca$hiv_status <- map.v1v3$hiv_status # cohort to split frames
pca$tooth_health <- map.v1v3$tooth_health # treatment
pca$visit_num <- map.v1v3$visit_num # time point data

# compute volatility -- euclidean distance over CLR transformed count data (aitchison distance between pairs)
# first need to melt our distance matrix to get pairs of samples and their distance
library(reshape2)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE", "FALSE")
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
library(dplyr)
group_by(dist_diff, V2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
)

# are they different?
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HEU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HEU",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
```

What about in just healthy teeth in diseased mouths?

```R






# filter samples

ps.dat.hcf <- subset_samples(ps.dat, aliquot_type=="H-CF")
# pull metadata table
map <- sample_data(ps.dat.hcf)
# this should look similar to the processing steps above. Basically want to create a new ASV table with duplicate samples across visit one and visit three
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
# pca 
counts.pca <- prcomp(counts.dist)
# extract the amount of variance for first four components
pc1 <- round(counts.pca$sdev[1]^2/sum(counts.pca$sdev^2),4) * 100
pc2 <- round(counts.pca$sdev[2]^2/sum(counts.pca$sdev^2),4) * 100
pc3 <- round(counts.pca$sdev[3]^2/sum(counts.pca$sdev^2),4) * 100
pc4 <- round(counts.pca$sdev[4]^2/sum(counts.pca$sdev^2),4) * 100

# extract scores for each sample for first four components
pca <- data.frame(PC1 = counts.pca$x[,1], 
                  PC2 = counts.pca$x[,2], 
                  PC3 = counts.pca$x[,3], 
                  PC4 = counts.pca$x[,4])

# add metadata to the pca scores
pca$studyID_FDI <- map.v1v3$studyID_FDI # individual tooth identifier
pca$hiv_status <- map.v1v3$hiv_status # cohort to split frames
pca$tooth_health <- map.v1v3$tooth_health # treatment
pca$visit_num <- map.v1v3$visit_num # time point data

# compute volatility -- euclidean distance over CLR transformed count data (aitchison distance between pairs)
# first need to melt our distance matrix to get pairs of samples and their distance
library(reshape2)
counts.melt <- melt(as.matrix(counts.dist))
# modify sample names in columns
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
# clean up data
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE", "FALSE")
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
library(dplyr)
group_by(dist_diff, V2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
)

# are they different?
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HI",]$value, dist_diff[dist_diff$V2 == "HEU",]$value)
wilcox.test(dist_diff[dist_diff$V2 == "HEU",]$value, dist_diff[dist_diff$V2 == "HUU",]$value)
```




QUESTION STILL OUTSTANDING -- ARE DIFFERENT HIV GROUPS MORE OR LESS LIKELY TO HAVE TOOTH HEALTH CHANGE? HOW STABLE IS HEALTHY TOOTH VS DISEASED TOOTH (WITHIN AND BETWEEN ALL SAMPLE GROUPS)
