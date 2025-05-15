# barcharts for volitility figure

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
map <- sample_data(ps.dat)
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
counts <- data.frame(table(map.v1v3$studyID_FDI))
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq > 1],]
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
counts <- otu_table(ps.dat)[row.names(map.v1v3),]
counts <- counts[,colSums(counts) > 0]
counts <- apply(counts,c(1,2),function(x) as.numeric(as.character(x)))
counts.clr <- clr(counts)
counts.dist <- dist(counts.clr, method="euclidean")
counts.melt <- melt(as.matrix(counts.dist))
counts.melt$Var1 <- gsub("V.*PQ", ".", counts.melt$Var1)
counts.melt$Var2 <- gsub("V.*PQ", ".", counts.melt$Var2)
counts.melt <- counts.melt[counts.melt$value != 0, ]
counts.melt$match <- ifelse(counts.melt$Var1 == counts.melt$Var2, "TRUE","FALSE")
matches <- counts.melt[counts.melt$match == "TRUE",]
matches <- as.data.frame(matches[,-2])
matches <- matches[!duplicated(matches),]
temp <- as.data.frame(cbind(map.v1v3$studyID_FDI, map.v1v3$hiv_status))
temp <- temp[!duplicated(temp),]
dist_diff <- merge(x=matches, y=temp, by.x="Var1", by.y="V1")
group_by(dist_diff, V2) %>%
summarise(
count = n(),
median = median(value, na.rm = TRUE),
SD = sd(value, na.rm = TRUE)
)


### HI
temp <- as_tibble(map.v1v3[map.v1v3$hiv_status == "HI",])
temp <- temp %>%
select(studyID_FDI, visit_num, tooth_health) %>%
pivot_wider(names_from=visit_num, values_from=tooth_health)
colnames(temp) <- c("studyID_FDI", "visit1", "visit3")
temp$change <- paste(temp$visit1, "to", temp$visit3, sep="_")
temp <- na.omit(temp)
temp <- temp %>% filter_all(all_vars(. !="NULL"))
changelist <- merge(x=dist_diff, y=temp, by.x="Var1", by.y="studyID_FDI")
changelist <- changelist %>%
mutate(change = recode(change,
'D_to_D' = 'No change disease',
'D_to_H' = 'Changed health state',
'E_to_D' = 'Changed health state',
'E_to_E' = 'No change disease',
'E_to_H' = 'Changed health state',
'H_to_D' = 'Changed health state',
'H_to_E' = 'Changed health state',
'H_to_H' = 'No change healthy'))
changelist$change <- factor(changelist$change, levels=c("Changed health state","No change disease", "No change healthy"))

library(ggplot2)
library(dplyr)

df <- changelist
avg_values <- df %>%
group_by(change) %>%
summarise(avg_value = mean(value))
global_avg_value <- mean(df$value)
total_count <- nrow(df)
options(repr.plot.width = 8, repr.plot.height =5)


p <- ggplot(df, aes(x = change)) +
# Bar chart for proportion of 'change'
geom_bar(aes(y = ..count.. / total_count), stat = "count", fill = "#8213A0") +
# "HI"="#8213A0", "HEU"="#FA78FA", "HUU"="#40A0FA"
# Line graph for scaled average 'value'
geom_line(data = avg_values, aes(x = change, y = avg_value / max(avg_value), group = 1),
color = "red", linewidth = 1) +
geom_point(data = avg_values, aes(x = change, y = avg_value / max(avg_value)),
color = "red", size = 2) +
# Add line for global average
geom_hline(yintercept = global_avg_value / max(avg_values$avg_value),
linetype = "dashed", color = "blue", size = 1) +
# Add labels with boxes to average points, offset so they don't cover the points
geom_label(data = avg_values,
aes(x = change,
y = avg_value / max(avg_values$avg_value),
label = round(avg_value, 1)), # Round to 1 decimal place
nudge_y = 0.05, # Offset vertically
label.size = 0.5, # Thickness of the box
fill = "white", # Background color of the box
color = "black") + # Color of the text
# Labels and title
labs(y = "Proportion / Scaled Average Value",
title = "Proportion of 'change' with Scaled Average 'value' and Global Average") +
ylim(0, 1) +
scale_y_continuous(
sec.axis = sec_axis(~ . * 25, name = "Average Value", breaks = seq(0, 25, by = 5)),
name = "Proportion"
) +
coord_flip() +
theme_minimal()
p
pdf("newvol_plot.HI.pdf")
p
dev.off()


### HEU
temp <- as_tibble(map.v1v3[map.v1v3$hiv_status == "HEU",])
temp <- temp %>%
select(studyID_FDI, visit_num, tooth_health) %>%
pivot_wider(names_from=visit_num, values_from=tooth_health)
colnames(temp) <- c("studyID_FDI", "visit1", "visit3")
temp$change <- paste(temp$visit1, "to", temp$visit3, sep="_")
temp <- na.omit(temp)
temp <- temp %>% filter_all(all_vars(. !="NULL"))
changelist <- merge(x=dist_diff, y=temp, by.x="Var1", by.y="studyID_FDI")
changelist <- changelist %>%
mutate(change = recode(change,
'D_to_D' = 'No change disease',
'D_to_H' = 'Changed health state',
'E_to_D' = 'Changed health state',
'E_to_E' = 'No change disease',
'E_to_H' = 'Changed health state',
'H_to_D' = 'Changed health state',
'H_to_E' = 'Changed health state',
'H_to_H' = 'No change healthy'))
changelist$change <- factor(changelist$change, levels=c("Changed health state","No change disease", "No change healthy"))

df <- changelist
avg_values <- df %>%
group_by(change) %>%
summarise(avg_value = mean(value))
global_avg_value <- mean(df$value)
total_count <- nrow(df)
options(repr.plot.width = 8, repr.plot.height =5)


p <- ggplot(df, aes(x = change)) +
# Bar chart for proportion of 'change'
geom_bar(aes(y = ..count.. / total_count), stat = "count", fill = "#FA78FA") +
# "HI"="#8213A0", "HEU"="#FA78FA", "HUU"="#40A0FA"
# Line graph for scaled average 'value'
geom_line(data = avg_values, aes(x = change, y = avg_value / max(avg_value), group = 1),
color = "red", linewidth = 1) +
geom_point(data = avg_values, aes(x = change, y = avg_value / max(avg_value)),
color = "red", size = 2) +
# Add line for global average
geom_hline(yintercept = global_avg_value / max(avg_values$avg_value),
linetype = "dashed", color = "blue", size = 1) +
# Add labels with boxes to average points, offset so they don't cover the points
geom_label(data = avg_values,
aes(x = change,
y = avg_value / max(avg_values$avg_value),
label = round(avg_value, 1)), # Round to 1 decimal place
nudge_y = 0.05, # Offset vertically
label.size = 0.5, # Thickness of the box
fill = "white", # Background color of the box
color = "black") + # Color of the text
# Labels and title
labs(y = "Proportion / Scaled Average Value",
title = "Proportion of 'change' with Scaled Average 'value' and Global Average") +
ylim(0, 1) +
scale_y_continuous(
sec.axis = sec_axis(~ . * 25, name = "Average Value", breaks = seq(0, 25, by = 5)),
name = "Proportion"
) +
coord_flip() +
theme_minimal()
p
pdf("newvol_plot.HEU.pdf")
p
dev.off()