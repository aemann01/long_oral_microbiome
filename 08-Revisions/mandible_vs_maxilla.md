```R
library(phyloseq, verbose=F)
library(ggplot2, verbose=F)
library(compositions, verbose=F)
library(reshape2, verbose=F)
library(tidyr, verbose=F)
library(dplyr, verbose=F)
library(ggdist, verbose=F)
library(ggpubr, verbose=F)
library(vegan)

load("../02-diversity/master_phyloseq.RData")
set.seed(48793)
sample_data(ps.dat)$tooth_type <- sample_data(ps.dat)$tooth_type %>% replace_na("unknown")
ps.dat <- subset_samples(ps.dat, tooth_type!="unknown")
ps.dat <- subset_samples(ps.dat, tooth_location!="unknown")
ps.dat <- subset_samples(ps.dat, tooth_age=="adult")
sample_data(ps.dat)$tooth_type[sample_data(ps.dat)$tooth_type == "premolar"] <- "molar"
ps.dat <- subset_samples(ps.dat, tooth_type == "molar")
ps.dat <- subset_samples(ps.dat, aliquot_type=="H-CF")
ps.dat
# convert to just lower and upper
sample_data(ps.dat)$tooth_location[sample_data(ps.dat)$tooth_location == "lower_left"] <- "lower"
sample_data(ps.dat)$tooth_location[sample_data(ps.dat)$tooth_location == "lower_right"] <- "lower"
sample_data(ps.dat)$tooth_location[sample_data(ps.dat)$tooth_location == "upper_left"] <- "upper"
sample_data(ps.dat)$tooth_location[sample_data(ps.dat)$tooth_location == "upper_right"] <- "upper"
table(sample_data(ps.dat)$tooth_location)

# first try with all samples
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
# are there signicant differences in the mandible vs maxilla samples?
map <- as(sample_data(ps.dat.clr), "data.frame")
adonis2(distance(ps.dat.clr, method="euclidean") ~ tooth_location, data = map)
# adonis2(formula = distance(ps.dat.clr, method = "euclidean") ~ tooth_location, data = map)
#           Df SumOfSqs      R2      F Pr(>F)
# Model      1     7598 0.00169 1.2316  0.002 **
# Residual 728  4491452 0.99831
# Total    729  4499050 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


################what about subsampled and by HIV status group?

# Initialize vectors to store results
p_values <- numeric(100)
r2_values <- numeric(100)

# first for HI
split <- subset_samples(ps.dat.clr, hiv_status=="HI")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_location))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_location) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_location, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
# [1] 0.4516753
# [1] 0.03196051

#### HEU
# Initialize vectors to store results
p_values <- numeric(100)
r2_values <- numeric(100)

split <- subset_samples(ps.dat.clr, hiv_status=="HEU")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_location))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_location) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_location, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
# [1] 0.3619232
# [1] 0.03011512

#### HUU
# Initialize vectors to store results
p_values <- numeric(100)
r2_values <- numeric(100)

split <- subset_samples(ps.dat.clr, hiv_status=="HUU")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_location))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_location) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_location, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
# [1] 0.1464465
# [1] 0.02800236
