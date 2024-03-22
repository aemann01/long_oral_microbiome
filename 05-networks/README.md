### Microbial coassociation networks

First load up conda environment on hillary

```bash
# running on hillary because of memory errors on palmetto
cd /home/allie/long_oral_microbiome/05-networks
conda activate 2023-Long_oral_microbiome
```

Set up R environment 

```R 
# Required packages
# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("circlize")
# install.packages("intergraph")
# install.packages("GGally")
# install.packages("sna")
# Install NetCoMi
# devtools::install_github("stefpeschel/NetCoMi", 
#                          dependencies = c("Depends", "Imports", "LinkingTo"),
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))
# load packages
library(phyloseq, verbose=F)
library(NetCoMi, verbose=F)
library(tidyverse, verbose=F)
library(reshape2, verbose=F)
library(igraph, verbose=F)
library(microbiome, verbose=F)
library(circlize, verbose=F)
library(intergraph, verbose=F)
library(network, verbose=F)
library(GGally, verbose=F)
library(sna, verbose=F)
```

Load data, collapse at roughly species level

```R
# load data
load("../02-diversity/master_phyloseq.RData")
ps.dat
ps.dat <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
ps.dat
```

Change network building parameters

```R
minobs <- 10
minprev <- 0.01
```

Networks within HIV status at each visit 

```R
# to clean up this code, need to run in a loop (and remove individual network output)
vars <- c("HUU", "HEU", "HI")
visits <- c("1", "2", "3")

for(i in 1:length(vars)){
	for(j in 1:length(visits)){
		sub <- subset_samples(ps.dat, hiv_status == vars[i] & visit_num == visits[j])
		print(paste("Running network for ", vars[i], " at visit ", visits[j]))
		print(paste("Number of samples: ", as.character(dim(otu_table(sub))[1])))
		print(paste("Number of ASVs: ", as.character(dim(otu_table(sub))[2])))
		# run spiec.easi
		pargs1 <- list(rep.num = 100, ncores = 25, seed = 10010)
		spiec.out <- spiec.easi(list(sub), method = 'mb', nlambda = 55, lambda.min.ratio=1e-1, pulsar.params = pargs1)
		spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(sub)))
		# write spiec-easi graph to file for anuran
		write.table(igraph::as_data_frame(spiec.graph, what="edges"), file=paste("spieceasi.ncol.", vars[i], "_", visits[j], ".txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
}
```

Tooth health specific networks at each visit

```R
vars <- c("H", "E", "D")
visits <- c("1", "2", "3")

for(i in 1:length(vars)){
	for(j in 1:length(visits)){
		sub <- subset_samples(ps.dat, tooth_health == vars[i] & visit_num == visits[j])
		print(paste("Running network for ", vars[i], " at visit ", visits[j]))
		print(paste("Number of samples: ", as.character(dim(otu_table(sub))[1])))
		print(paste("Number of ASVs: ", as.character(dim(otu_table(sub))[2])))
		# run spiec.easi
		pargs1 <- list(rep.num = 100, ncores = 25, seed = 10010)
		spiec.out <- spiec.easi(list(sub), method = 'mb', nlambda = 55, lambda.min.ratio=1e-1, pulsar.params = pargs1)
		spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(sub)))
		# write spiec-easi graph to file for anuran
		write.table(igraph::as_data_frame(spiec.graph, what="edges"), file=paste("spieceasi.ncol.", vars[i], "_", visits[j], ".txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
}
```

This one is a little more complicated as I want networks for each HIV status group at a specific tooth health level across all three visits.

```R
vars <- c("HUU", "HEU", "HI")
vars2 <- c("H", "D") # not doing for E because of low sample size
visits <- c("1", "2", "3")

for(i in 1:length(vars)){
	for(j in 1:length(vars2)){
		for(k in 1:length(visits)){
		sub <- subset_samples(ps.dat, hiv_status == vars[i] & tooth_health == vars2[j] & visit_num == visits[k])
		print(paste("Running network for ", vars[i], " and ", vars2[j], " at visit ", visits[k]))
		print(paste("Number of samples: ", as.character(dim(otu_table(sub))[1])))
		print(paste("Number of ASVs: ", as.character(dim(otu_table(sub))[2])))
		# run spiec.easi
		pargs1 <- list(rep.num = 100, ncores = 25, seed = 10010)
		spiec.out <- spiec.easi(list(sub), method = 'mb', nlambda = 55, lambda.min.ratio=1e-1, pulsar.params = pargs1)
		spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(sub)))
		# write spiec-easi graph to file for anuran
		write.table(igraph::as_data_frame(spiec.graph, what="edges"), file=paste("spieceasi.ncol.", vars[i], "_", vars2[j],"_", visits[k], ".txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
			}
		}
}

# getting a memory mapped error in this loop at HUU D 3 -- try manual network generation with fewer cores
sub <- subset_samples(ps.dat, hiv_status == "HUU" & tooth_health == "D" & visit_num == "3")
print(paste("Number of samples: ", as.character(dim(otu_table(sub))[1])))
print(paste("Number of ASVs: ", as.character(dim(otu_table(sub))[2])))
# run spiec.easi
pargs1 <- list(rep.num = 100, ncores = 10, seed = 10010)
spiec.out <- spiec.easi(list(sub), method = 'mb', nlambda = 55, lambda.min.ratio=1e-1, pulsar.params = pargs1)
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(sub)))
# write spiec-easi graph to file for anuran
write.table(igraph::as_data_frame(spiec.graph, what="edges"), file=paste("spieceasi.ncol.", "HUU", "_", "D","_", "3", ".txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# ok so there were some warnings with this one, try to run the rest of the HIV groups

vars <- c("HEU", "HI")
vars2 <- c("H", "D") # not doing for E because of low sample size
visits <- c("1", "2", "3")

for(i in 1:length(vars)){
	for(j in 1:length(vars2)){
		for(k in 1:length(visits)){
		sub <- subset_samples(ps.dat, hiv_status == vars[i] & tooth_health == vars2[j] & visit_num == visits[k])
		print(paste("Running network for ", vars[i], " and ", vars2[j], " at visit ", visits[k]))
		print(paste("Number of samples: ", as.character(dim(otu_table(sub))[1])))
		print(paste("Number of ASVs: ", as.character(dim(otu_table(sub))[2])))
		# run spiec.easi
		pargs1 <- list(rep.num = 100, ncores = 25, seed = 10010)
		spiec.out <- spiec.easi(list(sub), method = 'mb', nlambda = 55, lambda.min.ratio=1e-1, pulsar.params = pargs1)
		spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(sub)))
		# write spiec-easi graph to file for anuran
		write.table(igraph::as_data_frame(spiec.graph, what="edges"), file=paste("spieceasi.ncol.", vars[i], "_", vars2[j],"_", visits[k], ".txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
			}
		}
}

```

Now from each of our different network levels, generate core association networks using anuran

```bash
conda deactivate
# conda config --add channels ramellose
# conda create -n anuran pandas=1.5.3 anuran=1.1.0 
conda activate anuran
# create network folders, move speiceasi results to appropriate folder
mkdir global_network net_huu net_heu net_hi net_H net_E net_D net_huu_H net_huu_E net_huu_D net_heu_H net_heu_E net_heu_D net_hi_H net_hi_E net_hi_D
# first for our global network
cp spieceasi.ncol.huu.v* spieceasi.ncol.hi.v* spieceasi.ncol.heu.v* global_network







 all_networks
mv spieceasi.ncol.huu.v* networks_huu
mv spieceasi.ncol.heu.v* networks_heu
mv spieceasi.ncol.hi.v* networks_hi
# associations need to be shared by all of the networks within a group
anuran -i networks_hi networks_heu networks_huu -o CAN -draw -perm 5 -nperm 10 -size 1 -c -compare TRUE
# global network from all can networks from each visit (association must be present in two of the three to be included)
anuran -i all_networks -o full_CAN -draw -perm 5 -nperm 10 -size 0.6 -c 
# convert graphml to edge list
grep "<edge " full_CAN_all_networks_0.6_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > full_CAN_all_networks_0.6_intersection.txt
```

Additional conda packages for network analysis and reload our original conda environment

```bash
conda deactivate
conda activate 2023-Long_oral_microbiome
# conda install conda-forge::r-concaveman
# conda install conda-forge::r-ggraph
```


