# control for age and sex in relevant statistics

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
ps.dat
# pull metadata table
map <- sample_data(ps.dat)
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")   
library(compositions)

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


###########now collapse at species level, save copy to save on time\
ps.dat.L7 <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
###########


glom <- ps.dat.L7
# Does the balance of taxa on HIV status depend on sample age? Use the modal age to test
# table(sample_data(glom)$age_y)
#  10  11   3   4   5   6   7   8   9
# 236  10  42 196 250 224 364 339 299
# as most children are seven, let's go with that
# glom <- subset_samples(glom, age_y=="7")
# remove any taxa with fewer than 500 counts and in at least 1% of samples post merging 
glom <- filter_taxa(glom, function(x) sum(x > 500) > (0.01*length(x)), TRUE)
# pull data
dat <- as.data.frame(otu_table(glom))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# first need to get comparable variable to merge our data on
# head(dist_diff)
map$Var1 <- paste(map$study_id, map$FDI_code, sep=".")
# head(map)
# head(dist_diff)
# merge dist_diff with map
map <- merge(as.data.frame(map), dist_diff, by="Var1")

# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names


# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# x <- datmerge[,2:163]

# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)


# install.packages("coda4microbiome")
library(coda4microbiome)
x <- as.matrix(x) # ASV abundance table
x_time <- as.numeric(datmerge$visit_num) # time point
subject_id <- datmerge$study_id # subject id
y <- datmerge$value # response value (here is volatility)
ini_time <- 1
end_time <- 3
set.seed(852)
bal <- coda_glmnet_longitudinal(x, y, x_time, subject_id, ini_time, end_time, nfolds=10, showPlots = FALSE)