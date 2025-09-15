# 1. Build bowtie2 db from rpoC from HOMD genome
```sh
cd ~/rna_dohmain/rpoc/database
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ./ALL_genomes.fnn
grep "DNA-directed RNA polymerase subunit beta'" -A 1 ./ALL_genomes.fnn > rpoC_homd.fa
bowtie2-build rpoC_homd.fa rpoC_homd.db
```
# 2. Map the samples to rpoC homd db
```sh
cd ~/long_oral/12-align
mkdir bowtie2
ls ~/long_oral/01-read_processing/raw/filtN/*R1*.fastq.gz | sed 's/_S.*1.fastq.gz//' | sed 's/.*DM/DM/' | grep -v PCR | while read line; do
    fq1=$(ls ~/long_oral/01-read_processing/raw/filtN/${line}*R1_001.fastq.gz)
    fq2=$(ls ~/long_oral/01-read_processing/raw/filtN/${line}*R2_001.fastq.gz)
    bowtie2 -x ~/rna_dohmain/rpoc/database/rpoC_homd.db \
    -1 "$fq1"\
    -2 "$fq2"\
    --end-to-end  \
    --qc-filter \
    --no-unal \
    -t \
    --no-mixed \
    --threads 60 \
    -I 0 \
    -X 600 \
    -S ./bowtie2/$line.sam \
    2> ./bowtie2/$line.out \
    1> ./bowtie2/$line.err;
    done

parallel -j 60 'samtools view -bh {} > bowtie2/{/.}.bam' ::: bowtie2/*.sam 
parallel -j 60 'samtools sort -o bowtie2/{/.}.sorted.bam {}' ::: bowtie2/*.bam 
parallel -j 60 'samtools index {}' ::: bowtie2/*.sorted.bam 
```
# 3. Run feature counts
```sh
mkdir bowtie2/featurecounts
cd bowtie2
ls *.sorted.bam | parallel -j 60 '
    sample_name=$(basename {} .sorted.bam)
    out=./featurecounts/${sample_name}.out
    echo -e "Geneid\t${sample_name}" > "$out"
    samtools idxstats {} | awk "{print \$1, \$3}" | awk "{OFS=\"\t\"; print \$1, \$2}" >> "$out"
'
cd featurecounts
ulimit -n 4096
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=2;i<=NF;i+=2) printf "%s\t", $i; printf "\n"}' > read_counts.txt

cd ~/rna_dohmain/rpoc/trees
```
# 4. Run deseq2 on rpoC data
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(treeio)
setwd("~/rna_dohmain/rpoc/trees")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("~/long_oral/12-align/bowtie2/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
genecounts$X <- NULL
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
# filter gene count
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
# sum togther ads genes
# Create a new column with the part before the underscore
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.05:  99"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 74, 5.1%
# LFC < 0 (down)     : 25, 1.7%
# outliers [1]       : 291, 20%
# low counts [2]     : 500, 34%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  91"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 158, 11%
# LFC < 0 (down)     : 33, 2.3%
# outliers [1]       : 291, 20%
# low counts [2]     : 416, 29%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
#get baseman
# basemean_norm=rowMeans(counts(se_star, normalized=TRUE))


combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$SEQ_ID
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
# filter gene count
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# health is positive, dentin cavity negative
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HUUvHI")
head(resdf)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$SEQ_ID <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$SEQ_ID<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(resdf), as.data.frame(resHEUvHI), by = "SEQ_ID")
dfHIV <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "SEQ_ID")


#now do it for tooth health
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
# filter gene count
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
# sum togther ads genes
# Create a new column with the part before the underscore
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "E", "H"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.05:  99"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 74, 5.1%
# LFC < 0 (down)     : 25, 1.7%
# outliers [1]       : 291, 20%
# low counts [2]     : 500, 34%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  91"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 158, 11%
# LFC < 0 (down)     : 33, 2.3%
# outliers [1]       : 291, 20%
# low counts [2]     : 416, 29%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$SEQ_ID
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resEvD <- lfcShrink(se_star, coef="tooth_health_E_vs_D", type="apeglm")


# get heu v huu
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E",]
# filter gene count
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c( "E", "H"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# health is positive, dentin cavity negative
resHvE <- lfcShrink(se_star, coef="tooth_health_H_vs_E", type="apeglm")

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HvD")
head(resdf)

colnames(resEvD) <- paste0(colnames(resEvD), "_EvD")
resEvD$SEQ_ID <- rownames(resEvD)
head(resEvD)

colnames(resHvE) <- paste0(colnames(resHvE), "_HvE")
resHvE$SEQ_ID<- rownames(resHvE)
head(resHvE)

combined <- left_join(as.data.frame(resdf), as.data.frame(resEvD), by = "SEQ_ID")
dftooth <- left_join(as.data.frame(combined), as.data.frame(resHvE), by = "SEQ_ID")

dfHIV2 <- dfHIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
dftooth2 <- dftooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')
joined_df <- left_join(dfHIV, dftooth, by = c("SEQ_ID","Species", "genus", "Genus", "Genus_Species"))
# joined_df <- left_join(dfHIV, dftooth, by = c("SEQ_ID","Species", "genus", "Genus"))
# joined_df$Genus_Species <- paste0("Genus","_","Species")
write.table(joined_df, file="deseq_results_ADS-all_samples.txt", quote=F, sep="\t")
save.image("tree_data.all_samples.RData")
```
# 5. Uaw deseq2 to make tree
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)
setwd("/home/suzanne/rna_dohmain/rpoc/trees")
load("tree_data.all_samples.RData")

#make tree
tree <- read.tree("~/rna_dohmain/07-ads_expression/vince_trees/panaroo_0.6/core_genes/strep.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus == "Streptococcus")
filt_HIV $Genus_Species <- paste0(filt_HIV$Genus,"_",filt_HIV$Species)
filt_tooth <- filter(dftooth, Genus == "Streptococcus")
filt_tooth$Genus_Species <- paste0(filt_tooth$Genus,"_",filt_tooth$Species)
filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID', 'Genus_Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID', 'Genus_Species','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("SEQ_ID", "Genus_Species"))

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID, Genus_Species) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)
group_SEQ <- group_SEQ %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID)                  
  )

increased_by <- .005
first_tip <- 1.3
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(tree.root) %<+% group_SEQ +
  geom_tiplab(aes(label=species2), align = TRUE, size = 3.5, offset = 0.0) +
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.01,
         offset = 0.069,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 3,
         color = "black")

pdf("tree.strep_rpoC.all.pdf", height = 135, width =80)
  gheatmap(p1,
         group_SEQ3,
         width = 0.01,
         offset = 0.0825,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 3,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.strep_rpoC.all.pdf")


p1 <- gheatmap(p,
         log10(group_SEQ2),
         width = 0.01,
         offset = 0.069,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 3,
         color = "black")

pdf("tree.strep_rpoC.log10.all.pdf", height = 135, width =80)
  gheatmap(p1,
         log10(group_SEQ3),
         width = 0.01,
         offset = 0.0825,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 3,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.strep_rpoC.log10.all.pdf")
```
Look at distro of super star compared to other oralis
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/rpoc/")
homd <- read.csv("~/rna_dohmain/homd_map/annotations.updated.txt", header=T, sep="\t", quote="")
metadata <- read.table("~/long_oral/map_domhain_long_2.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("~/long_oral/12-align/bowtie2/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
genecounts$X <- NULL
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
# filter gene count
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
subcount <- subcount[, submap$sample_id]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
# sum togther ads genes
# Create a new column with the part before the underscore
#subset to just rpoC genes
rpoC_locus_tags <- homd$locus_tag[homd$gene_base == "rpoC"]
filtered_subcount <- subcount[rownames(subcount) %in% rpoC_locus_tags, ]

filtered_subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(filtered_subcount))
summed_subcount <- filtered_subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
#get relative abundance of dna
load("~/long_oral/master_phyloseq.RData")
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Streptococcus_oralis")
data <- psmelt(actino) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, visit_num, tooth_health, Abundance, V8, OTU)
red_dna$nucl <- "dna"
red_dna <- red_dna %>%
  rename(sample = Sample)
red_dna <- red_dna %>%
  rename(species = V8)
red_dna <- red_dna %>%
  rename(value = Abundance)
red_dna <- red_dna %>%
  rename(ASV = OTU)
red_dna[red_dna$sample == "DM00008V1PQ16-2", ]
#get relative abundance of rna
# red_rna <- as.data.frame(counts(se_star, normalized = TRUE)["SEQF6415.1", ])
# write results to file
strep_orals_ids <- unique((filter(homd, Corrected_Species == "Streptococcus_oralis")$SEQ_ID))
norm_counts <- as.data.frame(counts(se_star, normalized = FALSE))
strep_norm_counts_long <- norm_counts[rownames(norm_counts) %in% strep_orals_ids, ]
strep_norm_counts_long$SEQ_ID <- rownames(strep_norm_counts_long)
strep_norm_counts_long <- melt(strep_norm_counts_long, id.vars = "SEQ_ID", variable.name = "sample", value.name = "value")

red_rna <- strep_norm_counts_long %>%
  mutate(SEQ_ID = ifelse(SEQ_ID == "SEQF6415.1", SEQ_ID, "other_oralis")) %>%
  group_by(sample, SEQ_ID) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  arrange(sample, SEQ_ID)
red_rna$nucl <- "rna"

map$sample <- row.names(map)
meta <- as.data.frame(as.matrix(map)) %>% dplyr::select(sample, hiv_status, tooth_health, visit_num)
red_rna <- left_join(meta, red_rna, by = "sample")
red_rna$species <- "Streptococcus_oralis"
red_rna$ASV <- "NaN"
red_rna <- red_rna[, c("sample", "hiv_status", "visit_num", "tooth_health", "value", "species", "ASV" "nucl")]

#combine
combined_df <- rbind(red_rna, red_dna)
# plot
combined_df$hiv_status <- factor(combined_df$hiv_status, levels = c("HUU", "HEU", "HI"))
combined_df$tooth_health <- factor(combined_df$tooth_health, levels = c("H", "E", "D"))

rna_df <- combined_df[combined_df$nucl == "rna", ]
dna_df <- combined_df[combined_df$nucl == "dna", ]

rna_df <- red_rna %>%
  mutate(hiv_status = factor(hiv_status, levels = c("HUU", "HEU", "HI"))) %>%
  arrange(hiv_status, desc(value)) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

pdf("SEQF6415.DNA_rpoC.hiv.all.pdf", width =10, height =10)
ggplot(rna_df, aes(x = sample, y = value, fill = SEQ_ID)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1)
  ) +
  labs(
    title = "Streptococcus oralis RNA: SEQF6415.1 vs Other IDs (by HIV status) for mapped rpoC reads ",
    x = "Sample",
    y = "Raw Counts",
    fill = "SEQ_ID"
  ) +
  scale_y_continuous(labels = scales::comma)
dev.off()
system("~/.iterm2/imgcat ./SEQF6415.DNA_rpoC.hiv.all.pdf")

rna_df2 <- rna_df %>%
  mutate(tooth_health = factor(tooth_health, levels = c("H", "E", "D"))) %>%
  arrange(tooth_health, desc(value)) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

pdf("SEQF6415.DNA_rpoC.tooth.all.pdf", width =8, height =10)
ggplot(rna_df2, aes(x = sample, y = value, fill = SEQ_ID)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1)
  ) +
  labs(
    title = "Streptococcus oralis RNA: SEQF6415.1 vs Other IDs (by Tooth Health) for mapped rpoC reads",
    x = "Sample",
    y = "Raw Counts",
    fill = "SEQ_ID"
  ) +
  scale_y_continuous(labels = scales::comma)
dev.off()
system("~/.iterm2/imgcat ./SEQF6415.DNA_rpoC.tooth.all.pdf")
otu_data <- as.data.frame(as.matrix(otu_table(ps.dat)))
tax_data <- tax_table(ps.dat)
otu_df <- as.data.frame(t(as.data.frame(otu_data)))
tax_df <- as.data.frame(tax_data)

tax_df$FeatureID <- rownames(tax_df)
otu_df$FeatureID <- rownames(otu_df)
combined_data <- left_join(tax_df, otu_df, by = "FeatureID")
rownames(combined_data) <- combined_data$FeatureID
combined_data$FeatureID <- NULL
asv_numbers <- as.numeric(sub("ASV", "", rownames(combined_data)))
combined_data_ordered <- combined_data[order(asv_numbers), ]
head(combined_data_ordered)
write.table(combined_data_ordered, file = "ASVs_with_taxa.tsv", row.names = TRUE, sep = "\t", quote = FALSE)

```
# 7. Get unaligned read stats
```sh
cd ../
mkdir bowtie22
ls ~/rna_dohmain/rpoc/data/filtN/*R1*.fastq.gz | sed 's/_S.*1.fastq.gz//' | sed 's/.*DM/DM/' | while read line; do
    fq1=$(ls ~/rna_dohmain/rpoc/data/filtN/${line}*R1_001.fastq.gz)
    fq2=$(ls ~/rna_dohmain/rpoc/data/filtN/${line}*R2_001.fastq.gz)
    bowtie2 -x ./database/rpoC_homd.db \
    -1 "$fq1"\
    -2 "$fq2"\
    --end-to-end  \
    --qc-filter \
    -t \
    --no-mixed \
    --threads 60 \
    -I 0 \
    -X 600 \
    -S ./bowtie22/$line.sam \
    2> ./bowtie22/$line.out \
    1> ./bowtie22/$line.err;
    done

paste <(grep -nm1 "overall alignment" *out | awk '{print $1}' | sed 's/.out:13:/\t/g' | sed 's/%//') <(grep -nm1 "aligned concordantly 0" *out | awk '{print $3}' | sed 's/(//g' | sed 's/)//g' |sed 's/%//') <(grep -nm1 "aligned concordantly exactly" *out | awk '{print $3}' | sed 's/(//g' | sed 's/)//g' |sed 's/%//') <(grep -nm1 "aligned concordantly >1" *out | awk '{print $3}' | sed 's/(//g' | sed 's/)//g' |sed 's/%//') > alignment_stats.txt
echo "sample overall_alignment aligned_concordantly_0 aligned_concordantly_exactly_1 aligned_concordantly_>1" | cat - alignment_stats.txt > temp && mv temp alignment_stats.txt
sed -i 's/ /\t/g' alignment_stats.txt
```
# 8. Make alignment of strep reads and compare
```sh
awk '$8 == "Streptococcus_oralis" && $6 == "rpoC" { print $1 }' ~/rna_dohmain/homd_map/annotations.updated.txt | while read line; do grep -A 1 $line rpoC_homd.fa; done > strep_oralis.fna
mkdir vsearch && cd vsearch
cp ../strep_oralis.fna .
vsearch --cluster_fast strep_oralis.fna --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
# test using in silico pcr
sed 's/ /_/g' strep_oralis.fna > strep_oralis.fasta
perl in_silico_PCR.pl -s strep_oralis.fasta -a MAYGARAARMGNATGYTNCARGA -b GMCATYTGRTCNCCRTCRAA -f strep_oralis_pcr.fasta > results.txt
```
# 9. Make rpoC tree with ASVs for oralis
```sh
# get rpoC sequences
cd ~/rna_dohmain/rpoc
mkdir trees
cd trees
mkdir oralis_amp
cd oralis_amp
perl in_silico_PCR.pl -s ~/rna_dohmain/rpoc/database/vsearch/strep_oralis.fasta -a MAYGARAARMGNATGYTNCARGA -b GMCATYTGRTCNCCRTCRAA -e -f strep_oralis_pcr.fasta > results.txt
awk 'BEGIN { FS="\t" } NR==FNR && NR>1 { map[$1]=$2; next } /^>/ { sub(/^>/, ""); print ">" map[$1] } !/^>/' results.txt strep_oralis_pcr.fasta > renamed_pcr.fasta
grep -v '^amp_' renamed_pcr.fasta | sed 's/_.*//' > renamed_pcr.fna
# get ASV sequences
grep -w Streptococcus_oralis ../taxonomy_bac.txt | awk '{print $1}' |  while read line; do grep -wA 1 $line ../rep_set.fa; done > ASV.fna
# combine rpoC sequences and align
cat renamed_pcr.fna ASV.fna > all_oralis.fna
mafft --thread -1 all_oralis.fna > all_oralis.align.fna
# build tree
fasttree -nt all_oralis.align.fna > rpoc2.fast.tre
iqtree -s all_oralis.align.fna -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix oralis.rpoC -safe
```
Make tree in R now
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)
library(geiger)
setwd("/home/suzanne/rna_dohmain/rpoc/trees")
load("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/tree_data.RData")
homd <- read.csv("~/rna_dohmain/homd_map/annotations.updated.txt", header=T, sep="\t", quote="")
#make tree
tree <- read.tree("~/rna_dohmain/rpoc/trees/oralis_amp/oralis.rpoC.treefile")
# tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tree.root$tip.label <- gsub("^_R_", "", tree.root$tip.label)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# get the tip labels that have ADS for Strep
filt_HIV <- filter(dfHIV, Genus == "Streptococcus")
group_SEQ <- filt_HIV %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID)                  
  )

# make annotation file 
resdf <- as.data.frame(tree.root$tip.label)
names(resdf)[1] <- "SEQ_ID"
row.names(resdf) <- resdf$SEQ_ID
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
asv_ids <- resdf$SEQ_ID[grepl("^ASV", resdf$SEQ_ID)]
asv_ann <- data.frame(
  SEQ_ID = asv_ids,
  Genus = NA,
  Species = NA,
  Genus_Species = "ASV",
  Corrected_Species = "ASV"
)
rownames(asv_ann) <- asv_ann$SEQ_ID
ann2 <- rbind(ann, asv_ann)

rownames(ann2) <- ann$SEQ_ID
sortrow <- rownames(ann2)[order(match(rownames(resdf), rownames(ann2)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann2 <- ann2[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann2)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann2)
resdf <- resdf[ , -1]
resdf <- resdf %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID),
    species2 = str_remove(species2, "subsp\\.?_oralis_?"))
resdf$SEQ_ID

# check to see which have ADS
resdf <- resdf %>%
  mutate(has_ADS = SEQ_ID %in% filt_HIV$SEQ_ID)
resdf$type_indicator <- "non-ADS"
resdf$type_indicator[grepl("^ASV", resdf$SEQ_ID)] <- "ASV"
resdf$type_indicator[resdf$has_ADS == TRUE & !grepl("^ASV", resdf$SEQ_ID)] <- "ADS"

# add boostrapping
bs_values <- as.numeric(tree.root$node.label)
p <- ggtree(tree.root) %<+% resdf +
  geom_tiplab(aes(label=species2, color = type_indicator), align = TRUE, size = 3, offset = 0.0) +
  scale_color_manual(values = c("ADS" = "#F8766D", "non-ADS" = "gray", "ASV" ="#3CCCD0"))

pdf("tree.oralis_amp.pdf", height = 80, width = 60)
ggtree(tree.root) %<+% resdf +
  geom_tiplab(aes(label = species2, color = type_indicator), align = FALSE, size = 2, offset = 0.0) +
  scale_color_manual(values = c("ADS" = "#F8766D", "non-ADS" = "gray", "ASV" ="#3CCCD0"))+
  geom_text2(aes(subset = !isTip & !is.na(as.numeric(label)), label = label), hjust = -0.3, size = 3)
dev.off()
system("~/.iterm2/imgcat tree.oralis_amp.pdf")

# prune tree
tips_to_remove <- c(
  "ASV1489", "ASV3663", "ASV3563", "ASV4496", "ASV7039", "ASV768", "ASV1511",
  "ASV7846", "ASV7695", "ASV118", "ASV6911", "ASV5537", "ASV5597", "ASV1078",
  "ASV200", "ASV6016", "ASV1395", "ASV4259", "ASV1203", "ASV2374", "ASV343",
  "ASV1426", "ASV5298", "ASV6838", "ASV699", "ASV3823", "ASV7549", "ASV809",
  "ASV5821", "ASV1408", "ASV5481", "ASV2498", "ASV4734", "ASV4380", "ASV2592",
  "ASV3261", "ASV6728", "ASV3382", "ASV1698", "ASV7208", "ASV460", "ASV7781",
  "ASV6169", "ASV779", "ASV964", "ASV840", "ASV6163", "ASV418", "ASV4336",
  "ASV4536", "ASV191", "ASV2708", "ASV4879", "ASV7972", "ASV3014", "ASV4586",
  "ASV880", "ASV8411", "ASV5993", "ASV7978", "ASV889", "ASV7647", "ASV4162",
  "ASV6354", "ASV563", "ASV340", "ASV559", "ASV601", "ASV7766", "ASV3206",
  "ASV5768", "ASV4199", "ASV6836", "ASV6878", "ASV5271", "ASV6656", "ASV3144",
  "ASV3462"
)
tree_pruned <- drop.tip(tree.root, tips_to_remove)
tree.root2 <- midpoint_root(tree_pruned)

pdf("tree.oralis_amp.pruned.pdf", height = 40, width = 30)
ggtree(tree.root2) %<+% resdf +
  geom_tiplab(aes(label = species2, color = type_indicator), align = FALSE, size = 6, offset = 0.0) +
  scale_color_manual(values = c("ADS" = "#F8766D", "non-ADS" = "gray", "ASV" ="#3CCCD0"))+
  geom_text2(aes(subset = !isTip & !is.na(as.numeric(label)), label = label), hjust = -0.3, size = 3)
dev.off()
system("~/.iterm2/imgcat tree.oralis_amp.pruned.pdf")

```
# 10. Make large rpoC tree
```sh
cd ~/rna_dohmain/rpoc/trees
mkdir strep_rpoC
cd strep_rpoC
awk '$4 == "Streptococcus" && $6 == "rpoC" { print $1 }' ~/rna_dohmain/homd_map/annotations.updated.txt | while read line; do grep -A 1 $line ../../database/rpoC_homd.fa; done > strep.fna
vsearch --cluster_fast strep.fna --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
mafft --thread -1 strep.fna > strep.align.fna
iqtree -s strep.align.fna -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep_rpoC.core -safe
```

