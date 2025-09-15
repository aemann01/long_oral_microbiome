# 1. Run at ASV level
```R
library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$hiv_status, method="euclidean", padj = "fdr") # check signficane
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$aliquot_type, method="euclidean", padj = "fdr") # check signficane
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$tooth_health, method="euclidean", padj = "fdr") # check signficane

ps.dat.sub <- subset_samples(ps.dat, tooth_health == "D" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

save.image("HvD_amplicon.RData")
pdf("./diffab.HvD.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HvD.pdf")

#H HUUvHI
ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HUUvHI.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHI.asv.pdf")
ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")
ps.dat.sub2 <- subset_samples(ps.dat.sub2, tooth_health == "H" )
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HUUvHI.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHI.asv.pdf")
```
# 2. Run at species level
```R
library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])

ps.dat.sub <- subset_samples(glom, tooth_health == "D" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

save.image("HvD_amplicon.RData")
pdf("./diffab.HvD.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HvD.species.pdf")


ps.dat.sub <- subset_samples(glom, tooth_health == "D" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.HvD.species.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HvD.species.pdf")


ps.dat.sub <- subset_samples(glom, tooth_health == "E" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.EvH.species.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.EvH.species.pdf")

ps.dat.sub <- subset_samples(glom, tooth_health == "E" | tooth_health == "D")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.EvD.species.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.EvD.species.pdf")

#H HUUvHI
ps.dat.sub <- subset_samples(glom, hiv_status == "HUU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HUUvHI.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHI.species.pdf")
save.image("corn_species.RData")

ps.dat.sub <- subset_samples(glom, hiv_status == "HEU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HEUvHI.species.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HEUvHI.species.pdf")

ps.dat.sub <- subset_samples(glom, hiv_status == "HEU" | hiv_status == "HUU")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HEUvHUU.species.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HEUvHUU.species.pdf")

save.image("corn_species.RData")



# ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")
# ps.dat.sub2 <- subset_samples(ps.dat.sub2, tooth_health == "H" )
# da_analysis <- differentialTest(formula = ~ hiv_status,
#                                phi.formula = ~ hiv_status,
#                                formula_null = ~ 1,
#                                phi.formula_null = ~ hiv_status,
#                                test = "Wald",
#                                boot = FALSE,
#                                data = ps.dat.sub,
#                                fdr_cutoff = 0.05)
# pdf("./diffab.HUUvHI.asv.pdf", width = 20)
# plot(da_analysis, level=c("V8"))
# dev.off()
# system("~/.iterm2/imgcat ./diffab.HUUvHI.asv.pdf")
```
# 3. Run at ASV level
```R
library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/long_oral/corn_cob")
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])

ps.dat.sub <- subset_samples(ps.dat, tooth_health == "D" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.HvD.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HvD.asv.pdf")


ps.dat.sub <- subset_samples(ps.dat, tooth_health == "D" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.HvD.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HvD.asv.pdf")


ps.dat.sub <- subset_samples(ps.dat, tooth_health == "E" | tooth_health == "H")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.EvH.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.EvH.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, tooth_health == "E" | tooth_health == "D")
da_analysis <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)

pdf("./diffab.EvD.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.EvD.asv.pdf")

#H HUUvHI
ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HUUvHI.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHI.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HEU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HEUvHI.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HEUvHI.asv.pdf")

ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HEU" | hiv_status == "HUU")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
pdf("./diffab.HEUvHUU.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HEUvHUU.asv.pdf")

save.image("corn_ASV.RData")
```