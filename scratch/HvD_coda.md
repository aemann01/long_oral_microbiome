# 1. Species level
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat.sub <- subset_samples(glom, tooth_health == "D" | tooth_health == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.tooth.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.tooth.pdf")
```
H v E
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat.sub <- subset_samples(glom, tooth_health == "E" | tooth_health == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HvE.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HvE.pdf")
```
H v ED
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
sample_data(glom)$tooth_disease <- ifelse(sample_data(glom)$tooth_health == "H", "H", "ED")
ps.dat.sub <- subset_samples(glom, tooth_disease == "ED" | tooth_disease == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_disease)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HvED.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HvED.pdf")
```
D v E
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# sample_data(glom)$tooth_disease <- ifelse(sample_data(glom)$tooth_health == "H", "H", "ED")
ps.dat.sub <- subset_samples(glom, tooth_health == "E" | tooth_health == "D")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.EvD.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.EvD.pdf")
```
HI v HUU
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# sample_data(glom)$tooth_disease <- ifelse(sample_data(glom)$tooth_health == "H", "H", "ED")
ps.dat.sub <- subset_samples(glom, hiv_status == "HI" | hiv_status == "HUU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$hiv_status)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
# geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HIvHUU.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HIvHUU.pdf")

ps.dat.sub <- subset_samples(glom, hiv_status == "HI" | hiv_status == "HEU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# x# response variable
y <- as.factor(dat$hiv_status)
# y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
# geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HIvHEU.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HIvHEU.pdf")

ps.dat.sub <- subset_samples(glom, hiv_status == "HUU" | hiv_status == "HEU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# x# response variable
y <- as.factor(dat$hiv_status)
# y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
# geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HUUvHEU.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUUvHEU.pdf")
```
# 2. ASV level
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat.sub <- subset_samples(ps.dat, tooth_health == "D" | tooth_health == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.tooth.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.tooth.asv.pdf")
```
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat.sub <- subset_samples(ps.dat, tooth_health == "E" | tooth_health == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HvE.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HvE.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, tooth_health == "E" | tooth_health == "D")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.DvE.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.DvE.asv.pdf")


ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HI" | hiv_status == "HUU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$hiv_status)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HIvHUU.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HIvHUU.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HI" | hiv_status == "HEU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$hiv_status)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HIvHEU.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HIvHEU.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HEU")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$hiv_status)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HUUvHEU.asv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUUvHEU.asv.pdf")
```
H v ED
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
sample_data(ps.dat)$tooth_disease <- ifelse(sample_data(ps.dat)$tooth_health == "H", "H", "ED")
ps.dat.sub <- subset_samples(ps.dat, tooth_disease == "ED" | tooth_disease == "H")# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(ps.dat.sub , function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$tooth_disease)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HvED.asv.pdf", height =20)
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HvED.asv.pdf")
```