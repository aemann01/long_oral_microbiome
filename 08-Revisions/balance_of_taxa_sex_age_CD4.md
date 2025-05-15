```R
library(phyloseq, verbose=F)
library(ggplot2, verbose=F)
library(compositions, verbose=F)
library(reshape2, verbose=F)
library(tidyr, verbose=F)
library(dplyr, verbose=F)
library(ggdist, verbose=F)
library(ggpubr, verbose=F)

load("../02-diversity/master_phyloseq.RData")
ps.dat.L7 <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])

# Age
glom <- ps.dat.L7
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
dif <- dim(datmerge)[2] - dim(map)[2] + 1
x <- datmerge[,2:dif]
x <- as.matrix(x) # ASV abundance table

y <- as.integer(datmerge$age_y) # response value 

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_age.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_age.pdf")
bal$`predictions plot`
dev.off()


## Sex
glom <- ps.dat.L7
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
dif <- dim(datmerge)[2] - dim(map)[2] + 1
x <- datmerge[,2:dif]
x <- as.matrix(x) # ASV abundance table

y <- as.factor(datmerge$sex) # response value 

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_sex.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_sex.pdf")
bal$`predictions plot`
dev.off()

## CD4 counts
glom <- ps.dat.L7
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
dif <- dim(datmerge)[2] - dim(map)[2] + 1
x <- datmerge[,2:dif]
x <- as.matrix(x) # ASV abundance table

y <- as.factor(datmerge$sex) # response value 

set.seed(852)
bal <- coda_glmnet(x, y)

# save predictions plots
coef <- bal$`log-contrast coefficients`
positives <- which(coef>=0)
op <- order(coef[positives], decreasing = TRUE)
pdf("microbial_sig_sex.pdf")
bal$`signature plot`
dev.off()
pdf("microbial_pred_sex.pdf")
bal$`predictions plot`
dev.off()