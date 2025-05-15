# subsample biogeography analyses

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
set.seed(48793)
sample_data(ps.dat)$tooth_type <- sample_data(ps.dat)$tooth_type %>% replace_na("unknown")
ps.dat <- subset_samples(ps.dat, tooth_type!="unknown")
ps.dat <- subset_samples(ps.dat, tooth_age=="adult")
sample_data(ps.dat)$tooth_type[sample_data(ps.dat)$tooth_type == "premolar"] <- "molar"
sample_data(ps.dat)$tooth_type[sample_data(ps.dat)$tooth_type == "canine"] <- "incisor"
ps.dat

ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")

# HI samples first

set.seed(5933)
split <- subset_samples(temp, hiv_status=="HI")

# get number of samples per group
subnum <- min(table(sample_data(split)$tooth_type))
sample_df <- as.data.frame(sample_data(split))
subsampled_samples <- sample_df %>%
  group_by(tooth_type) %>%  # Replace `hiv_status` with your factor variable
  sample_n(min(subnum, n()))
df <- data.frame(subsampled_samples)
downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

# permanova
map <- as(sample_data(downsamp), "data.frame")
adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999

# adonis2(formula = distance(downsamp, method = "euclidean") ~ tooth_type, data = map)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1     4990 0.05852 0.9945  0.526
# Residual 16    80277 0.94148
# Total    17    85266 1.00000

# HUU
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
split <- subset_samples(temp, hiv_status=="HUU")
sample_df <- as.data.frame(sample_data(split))
# get number of samples per group
subnum <- min(table(sample_data(split)$tooth_type))

subsampled_samples <- sample_df %>%
  group_by(tooth_type) %>%  
  sample_n(min(subnum, n()))

df <- data.frame(subsampled_samples)
downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

# permanova
map <- as(sample_data(downsamp), "data.frame")
adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)
# adonis2(formula = distance(downsamp, method = "euclidean") ~ tooth_type, data = map)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1     8372 0.04958 1.3564  0.002 **
# Residual 26   160484 0.95042
# Total    27   168856 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# HEU
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
split <- subset_samples(temp, hiv_status=="HEU")
sample_df <- as.data.frame(sample_data(split))
# get number of samples per group
subnum <- min(table(sample_data(split)$tooth_type))

subsampled_samples <- sample_df %>%
  group_by(tooth_type) %>%  
  sample_n(min(subnum, n()))

df <- data.frame(subsampled_samples)
downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

# permanova
map <- as(sample_data(downsamp), "data.frame")
adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)
# adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999

# adonis2(formula = distance(downsamp, method = "euclidean") ~ tooth_type, data = map)
#          Df SumOfSqs     R2      F Pr(>F)
# Model     1     7564 0.0566 1.3199  0.001 ***
# Residual 22   126073 0.9434
# Total    23   133637 1.0000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Because these are variable -- let's do a permutation test

```R
# Initialize vectors to store results
p_values <- numeric(100)
r2_values <- numeric(100)

# first for HI
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
split <- subset_samples(temp, hiv_status=="HI")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_type))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_type) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
```

Permutation analysis for HUU

```R
p_values <- numeric(100)
r2_values <- numeric(100)

# first for HI
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
split <- subset_samples(temp, hiv_status=="HUU")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_type))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_type) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
```

For HEU

```R
p_values <- numeric(100)
r2_values <- numeric(100)

# first for HI
temp <- subset_samples(ps.dat.clr, aliquot_type=="H-CF")
split <- subset_samples(temp, hiv_status=="HEU")
sample_df <- as.data.frame(sample_data(split))
subnum <- min(table(sample_data(split)$tooth_type))

n_reps <- 100

for (i in 1:n_reps){
	subsampled_samples <- sample_df %>%
  		group_by(tooth_type) %>%  
  		sample_n(min(subnum, n()))

	df <- data.frame(subsampled_samples)
	downsamp <- subset_samples(split, file_R1 %in% df$file_R1)

	# permanova
	map <- as(sample_data(downsamp), "data.frame")
	permanova_result <- adonis2(distance(downsamp, method="euclidean") ~ tooth_type, data = map)

	p_values[i] <- permanova_result$`Pr(>F)`[1]
	r2_values[i] <- permanova_result$R2[1]
}

# adjust all p values using FDR
adjpval <- p.adjust(p_values, method = "fdr")

# calculate average p value and r2 value
mean(adjpval)
mean(r2_values)
```






