### Community fluctuations before and after high Strep mutans colonization

From previous data we've generated it doesn't look like there are much if any mutans on teeth in early stages of caries development, why is this? What communities are present before and after mutans takes over a tooth? Can we predict high strep mutans based on the community before it takes over? Are there predictable taxa that colonize after the community of mutans collapses (for whatever reason)?

Activate environment

```bash
# on hillary
cd /home/allie/long_oral_microbiome/06-mutans
conda activate 2023-Long_oral_microbiome
```

Install and load required libraries

```R
# install.packages("caTools")
# install.packages("ranger")
# install.packages("fastshap")
# install.packages("kernelshap")
# install.packages("iml")
# install.packages("vip")
# install.packages("shapviz")

library(phyloseq, verbose=F)
library(tidyverse, verbose=F)
library(caTools, verbose=F)
library(ranger, verbose=F)
library(fastshap, verbose=F)
library(kernelshap, verbose=F)
library(iml, verbose=F)
library(vip, verbose=F)
library(shapviz, verbose=F)
library(reshape2, verbose=F)
```

Load data and collapse roughly by species level

```R
load("../02-diversity/master_phyloseq.RData")
ps.dat
ps.dat <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat
```

First need to identify teeth with a high level of Strep mutans at any time point.

```R
# get the relative abundance of strep mutans for all samples at any time point
relabund <- transform_sample_counts(ps.dat, function(x) x / sum(x))
# which ASV ID is strep mutans?
tax <- as.data.frame(tax_table(relabund))
smut <- row.names(tax[str_detect(tax$V8, "Streptococcus_mutans"),])
abund <- as.data.frame(otu_table(relabund))
# pull relative abundance for Smutans data by sample ID
temp <- as.data.frame(cbind(row.names(abund), abund[, names(abund) %in% smut]))
colnames(temp) <- c("sample_id", "relabund_smutans")
# which samples have higher than 10% strep mutans?
highsmut <- temp[temp$relabund >= 0.1,]
# lower than 1% strep mutans
lowsmut <- temp[temp$relabund <= 0.05,]

# now need to create a new data frame with smut proportion, sample ID, study ID, FDI code, and visit number all separate for both low and high 
# individual ID
ind <- gsub("V.*", "", highsmut$sample_id)
# FDI, clean up
fdi <- gsub("D.*PQ", "", highsmut$sample_id)
fdi <- gsub("-.*$", "", fdi)
# visit 
visit <- gsub("D.*V", "", highsmut$sample_id)
visit <- gsub("PQ.*", "", visit)
# combine into single dataframe
highsmut <- cbind(highsmut, ind, fdi, visit)
# clean up headers so we don't get confused
colnames(highsmut) <- c("sample_id.high", "smut.high", "ind", "fdi", "visit.high")

# same for low smutans
# individual ID
ind <- gsub("V.*", "", lowsmut$sample_id)
# FDI, clean up
fdi <- gsub("D.*PQ", "", lowsmut$sample_id)
fdi <- gsub("-.*$", "", fdi)
# visit 
visit <- gsub("D.*V", "", lowsmut$sample_id)
visit <- gsub("PQ.*", "", visit)
# combine into single dataframe
lowsmut <- cbind(lowsmut, ind, fdi, visit)
# clean up headers so we don't get confused
colnames(lowsmut) <- c("sample_id.low", "smut.low", "ind", "fdi", "visit.low")

# merge high and low smut by individual ID and FDI code 
highlow.smut <- merge(highsmut, lowsmut, by.x=c("ind", "fdi"), by.y=c("ind", "fdi"))
# is low smutans sample before or after high smutans?
highlow.smut$before_after <- ifelse(highlow.smut$visit.high > highlow.smut$visit.low, "before", 
	ifelse(highlow.smut$visit.high < highlow.smut$visit.low, "after", "None"))

temp.df <- highlow.smut
# finally, add in high s mutans samples as "during"
temp.befaf <- as.data.frame(cbind(temp.df$sample_id.low, temp.df$ind, temp.df$fdi, temp.df$smut.low, temp.df$visit.low, as.character(temp.df$before_after)))
colnames(temp.befaf) <- c("sampleID", "ind", "fdi", "smut", "visit", "period")
# add in during samples
temp.dur <- as.data.frame(cbind(temp.df$sample_id.high, temp.df$ind, temp.df$fdi, temp.df$smut.high, temp.df$visit.high))
temp.dur <- cbind(temp.dur, rep("during", length(temp.dur$V1)))
colnames(temp.dur) <- c("sampleID", "ind", "fdi", "smut", "visit", "period")
temp.all <- rbind(temp.befaf, temp.dur)
# remove duplicate rows
temp.all <- temp.all[!duplicated(temp.all),]

## below code saving if this problem happens again (did with 10% vs 10% analysis)
# ## so one sample (at visit two) is both a before AND an after because smutans is high at visit one and visit three which mucks things up. this is a little unsatisfying but I'll include the one that is representative of the condition I have the least of
# table(temp.all$period) # there are fewer "befores" than any other group -- only keep the before entry
# n_occur <- data.frame(table(temp.all$sampleID))
# dups <- as.character(n_occur[n_occur$Freq > 1,][[1]])
# index <- row.names(temp.all[temp.all$sampleID %in% dups & temp.all$period == "before",])
# temp.all <- temp.all[-c(as.numeric(index)),]

# row names
row.names(temp.all) <- temp.all$sampleID
# data frame for random forest
df.rf <- temp.all

table(df.rf$period)

 # after before during
 #    26     19     36
```

Random forest analysis and Shapley Predictions

```R
set.seed(151)

# first need to format data for random forest predictions
asv_tab <- as.data.frame(otu_table(ps.dat)) # this is our ASV table after species level glom

df.rf$period <- factor(df.rf$period)
# filter to only include the samples we want to look at 
sub_asv <- asv_tab[row.names(asv_tab) %in% df.rf$sampleID,]
dim(sub_asv)
dim(df.rf)

# normalize and scale 
asv_tab_norm <- sweep(sub_asv, 2, colSums(sub_asv), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
# remove nas
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]

# set variable of interest
asv_tab_var$var <- df.rf[rownames(asv_tab_var), "period"]

######## create random forest model
fit <- ranger(var ~ ., data = asv_tab_var, scale.permutation.importance = TRUE, importance = 'permutation')
# confusion matrix
fit$confusion.matrix
#         predicted
# true     after before during
#   after     16      0     10
#   before     4      3     12
#   during     7      1     28
# variable importance
pdf("rf.period.importance.pdf")
vip(fit, title = "Variable Importance")
dev.off()
# ranger::importance(fit)
```

Post hoc explanatory variables for each mutans group

```R
# save.image("shap.RData")
# load("shap.RData")
# run in parallel
library(doParallel)
registerDoParallel(cores = 30)
set.seed(151)

# refit random forest for shapley explainations 
rfo <- ranger(var ~ ., data = asv_tab_var, probability = TRUE)

# first find predictive variables for during (sanity check, should be Strep mutans ASV11)
pfun <- function(object, newdata) {  # prediction wrapper
  predict(object, data = newdata)$predictions[,"during"]
}
X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.4439358
ex.during <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

shv.during <- shapviz(ex.during)
pdf("sv_importance.during.pdf")
sv_importance(shv.during)  
dev.off()

# works!!

# Predictive variables for before
pfun <- function(object, newdata) {  # prediction wrapper
  predict(object, data = newdata)$predictions[,"before"]
}
X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.238184
ex.before <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

# shapley based feature importance score -- mean of the absolute value of the feature contribution for each column 
shv.before <- shapviz(ex.before)
pdf("sv_importance.before.pdf")
sv_importance(shv.before)  
dev.off()








# Predictive variables for after
pfun <- function(object, newdata) {  # prediction wrapper
  predict(object, data = newdata)$predictions[,"after"]
}
X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.3174867
ex.after <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

# shapley based feature importance score -- mean of the absolute value of the feature contribution for each column 
shv.after <- shapviz(ex.after)
pdf("sv_importance.after.pdf")
sv_importance(shv.after)  
dev.off()
```



