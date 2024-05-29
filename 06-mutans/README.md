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
# we need to convert the relabund to numeric but don't make it scientific notation or everything will get all wierd later on
# turn off scientific notation
options(scipen = 999)
temp$relabund_smutans <- as.numeric(as.character(temp$relabund_smutans))

# which samples have higher than 10% strep mutans?
highsmut <- temp[temp$relabund >= 0.1,]
# lower than 5% strep mutans
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

df.rf$smutTEST <- as.numeric(as.character(df.rf$smut))

table(df.rf$period)

 # after before during
 #    15     12     24
```

Random forest analysis and Shapley Predictions

```R
set.seed(151)

# first need to format data for random forest predictions
asv_tab <- as.data.frame(otu_table(ps.dat)) # this is our ASV table after species level glom

df.rf$period <- factor(df.rf$period)
# filter to only include the samples we want to look at 
sub_asv <- asv_tab[row.names(asv_tab) %in% df.rf$sampleID,]
# also remove strep mutans as a variable -- this is the response for during -- let's see what happens when we remove it?
sub_asv <- sub_asv[, colnames(sub_asv) != "ASV11"]

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
#   after      7      0      8
#   before     2      4      6
#   during     3      1     20
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
  predict(object, data = newdata)$predictions[,"during"] # conditional probability that a sample is "during" high smut
}


X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.4692195 # baseline prediction that any sample is "during"
ex.during <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

shv.during <- shapviz(ex.during)
pdf("sv_importance.during.pdf")
sv_importance(shv.during)  
dev.off()

# bee plot
pdf("sv_importance.during.bee.pdf")
sv_importance(shv.during, kind="bee")
dev.off()

# works!!

# Predictive variables for before
pfun <- function(object, newdata) {  # prediction wrapper
  predict(object, data = newdata)$predictions[,"before"]
}
X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.2330286
ex.before <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

# shapley based feature importance score -- mean of the absolute value of the feature contribution for each column 
shv.before <- shapviz(ex.before)
pdf("sv_importance.before.pdf")
sv_importance(shv.before)  
dev.off()

# bee plot
pdf("sv_importance.before.bee.pdf")
sv_importance(shv.before, kind="bee")
dev.off()

# Predictive variables for after
pfun <- function(object, newdata) {  # prediction wrapper
  predict(object, data = newdata)$predictions[,"after"]
}
X <- subset(asv_tab_var, select = -var)
mean(pfun(rfo, newdata = asv_tab_var))
# [1] 0.2977519
ex.after <- explain(rfo, X = X, pred_wrapper = pfun, nsim = 100, adjust = TRUE,
                 shap_only = FALSE, parallel = TRUE)

# shapley based feature importance score -- mean of the absolute value of the feature contribution for each column 
shv.after <- shapviz(ex.after)
pdf("sv_importance.after.pdf")
sv_importance(shv.after)  
dev.off()

# bee plot
pdf("sv_importance.after.bee.pdf")
sv_importance(shv.after, kind="bee")
dev.off()
```

Line plot of smutans proportion by period

```R
df.rf$period <- as.factor(factor(df.rf$period, levels=c("before", "during", "after")))

df.summary <- df.rf %>%
  group_by(period) %>%
  summarise(
    sd = sd(smutTEST),
    smutTEST = mean(smutTEST)
  )
df.summary

pdf("smut_prop.pdf", width = 10, height = 3)
ggplot(df.rf, aes(period, smutTEST)) +
  geom_jitter(position = position_jitter(0.2)) + 
  geom_line(aes(group = 1), data = df.summary) +
  geom_errorbar(
    aes(ymin = smutTEST-sd, ymax = smutTEST+sd),
    data = df.summary, width = 0.2) +
  geom_point(data = df.summary, size = 3) +
  theme_minimal()
dev.off()
```

Using what we found in the shap analysis -- make the same figure but this time with top two taxa from after and before added in

```R
# filter samples from abundance table
sampsfilt <- abund[row.names(abund) %in% row.names(temp.all),]
# now filter ASVs
wantedIDs <- c("ASV11", "ASV54", "ASV298", "ASV478", "ASV141")
# pull relative abundance data for each ASV
temp <- as.data.frame(cbind(row.names(sampsfilt), sampsfilt[, names(sampsfilt) %in% wantedIDs]))
# join with temp.all with period info
merged <- merge(temp.all, temp, by='row.names')
row.names(merged) <- merged$Row.names
# leaving ASV11 in there as a sanity check, should match smut column
merged$period <- as.factor(factor(merged$period, levels=c("before", "during", "after")))
# melt dataframe 
melted <- melt(merged, id = c("period", "ind", "visit", "Row.names", "sampleID", "fdi", "smut", "row.names(sampsfilt)"))


df.summary <- melted %>%
  group_by(period,variable) %>%
  summarise(
    sd = sd(value),
    mean = mean(value)
  )
df.summary

pdf('test.pdf')
ggplot(df.summary, aes(period, mean)) +
  geom_col(aes(fill = variable), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = mean, ymax = mean+sd, group = period),
    width = 0.2, position = position_dodge(0.8)
    )
 dev.off()
```

For final bit of figure, make a cap scale plot showing total community composition of samples before, during, and after?

```R
# recreate our phyloseq object but with our strep mutans data as the sample data to subsample
# reload ASV level ps.dat object
load("../02-diversity/master_phyloseq.RData")
ps.dat
# subsample
sample_data(ps.dat) <- df.rf
ps.dat
# remove zero asvs
filt <- filter_taxa(ps.dat, function(x) sum(x) > 1, TRUE)
filt
# transform to even sampling depth 
filt <- transform_sample_counts(filt, function(x) 1E6 * x/sum(x))

ordcap <- ordinate(ps.dat, "PCoA", "unifrac", weighted=TRUE)
pdf("capscale_plt.period.pdf") 
plot_ordination(ps.dat, ordcap, "samples", color="period") + theme_minimal()
dev.off()
```

After looking at multiple ordination methods, there aren't any really clear clusters -- possibly because we have extreme outliers + these are the same teeth over time -- might come back to this

Now want to pull samples with high values for ASV54 (Strep sanguinis) before S. mutans -- is this a specific ASV or ASVs that are driving this? Should be a shap value of higher than 0.05

```R
before.vals <- get_shap_values(shv.before)
row.names(before.vals) <- row.names(get_feature_values(shv.before))
before.vals <- as.data.frame(before.vals)
# what are the samples with the highest shap values?
head(before.vals[order(before.vals$ASV54, decreasing=TRUE),])
# [1] "DM00024V1PQ54-1" "DM00325V1PQ84"   "DM00074V1PQ54"   "DM00393V2PQ61"
# [5] "DM00312V1PQ16"   "DM00092V1PQ55"
```

Now I want to see if these share specific ASVs

```R
# reload full uncollapsed data
load("../02-diversity/master_phyloseq.RData")
ps.dat
samps <- row.names(head(before.vals[order(before.vals$ASV54, decreasing=TRUE),]))
sub.ps <- prune_samples(sample_names(ps.dat) %in% samps, ps.dat)
# only focus on sanguinis ASVs
sub.ps <- subset_taxa(sub.ps, V8=="Streptococcus_sanguinis")
pdf("sub_barplot.pdf")
plot_bar(sub.ps, fill="V9")
dev.off()
```

These do not have very similar profiles -- not super informative

