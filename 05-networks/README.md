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

# Problems with D and E, re run with just H

vars <- c("HUU", "HEU", "HI")
vars2 <- c("H") # not doing for E OR D because of low sample size
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
mkdir global_network net_huu net_heu net_hi net_H net_E net_D net_huu_H net_heu_H  net_hi_H  
# first for our global network
cp spieceasi.ncol.HUU_1.txt spieceasi.ncol.HUU_2.txt spieceasi.ncol.HUU_3.txt spieceasi.ncol.HI_1.txt spieceasi.ncol.HI_2.txt spieceasi.ncol.HI_3.txt spieceasi.ncol.HEU_1.txt spieceasi.ncol.HEU_2.txt spieceasi.ncol.HEU_3.txt global_network
mv spieceasi.ncol.HUU_1.txt spieceasi.ncol.HUU_2.txt spieceasi.ncol.HUU_3.txt net_huu/
mv spieceasi.ncol.HEU_1.txt spieceasi.ncol.HEU_2.txt spieceasi.ncol.HEU_3.txt net_heu/
mv spieceasi.ncol.HI_1.txt spieceasi.ncol.HI_2.txt spieceasi.ncol.HI_3.txt net_hi/
mv spieceasi.ncol.D_* net_D
mv spieceasi.ncol.E_* net_E
mv spieceasi.ncol.H_* net_H
mv spieceasi.ncol.HEU_H_* net_heu_H
mv spieceasi.ncol.HUU_H_* net_huu_H
mv spieceasi.ncol.HI_H_* net_hi_H

# associations need to be shared by all of the networks within a group
anuran -i net_hi net_heu net_huu net_heu_H net_huu_H net_hi_H -o CAN -draw -perm 5 -nperm 10 -size 1 -c -compare TRUE
# global network from all can networks from each visit (association must be present in two of the three to be included)
anuran -i global_network -o full_CAN -draw -perm 5 -nperm 10 -size 0.6 -c 
# convert graphml to edge list
grep "<edge " full_CAN_global_network_0.6_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > full_CAN_global_network_0.6_intersection.txt
```

Additional conda packages for network analysis and reload our original conda environment

```bash
conda deactivate
conda activate 2023-Long_oral_microbiome
# conda install conda-forge::r-concaveman
# conda install conda-forge::r-ggraph
```

Detect community modularity from global CAN network

```R
set.seed(782469)
# install.packages("oaqc")
library(igraph)
library(ggforce)
library(concaveman)
library(ggraph)
library(oaqc)

dat <- read.table("full_CAN_all_networks_0.6_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
coords <- layout_nicely(can) 
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.7162675 # high modularity reflects dense connections within communities and sparse connections across communities
pdf("all_can_network.clust.labs.pdf")
plot(clust, can, layout=coords, vertex.size=7, vertex.color=V(can)$color)
dev.off()
# get dataframe of representative ASVs with their cluster membership
clusts <- list()
for(i in 1:length(clust[])){
     col1 <- eval(parse(text=paste0("clust[]$`",i,"`", sep="")))
     reps <- rep(i, length(col1))
     clusts[[i]] <- cbind(col1, reps)
     }
cdat <- as.data.frame(do.call(rbind, clusts))
colnames(cdat) <- c("repASV", "cluster")
# write clusters to file
write.table(cdat, file="clusters.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
```

Now want to populate our cluster dataframe with taxonomy and whether or not the ASVs are found in our HIV group specific CAN networks

```bash
awk '{print $1}' clusters.txt | grep -v "rep" | while read line; do grep -m 1 -w $line ../02-diversity/taxonomy_bac.filt.txt | awk -F"\t" '{print $9}'; done > temp
sed -i '1 i\taxonomy' temp
paste clusters.txt temp > temp2
# check if ASV ID exists in CAN network
awk '{print $1}' temp2 | grep -v "rep" | while read line; do if cat CAN_networks_huu_1_intersection.graphml | grep -w -m 1 -q $line; then
     echo "1"
else
     echo "0"
fi; 
done > HUU_presabs
sed -i '1 i\HUU_CAN' HUU_presabs
# do the same with the other two can networks
awk '{print $1}' temp2 | grep -v "rep" | while read line; do if cat CAN_networks_heu_1_intersection.graphml | grep -w -m 1 -q $line; then
     echo "1"
else
     echo "0"
fi; 
done > HEU_presabs
sed -i '1 i\HEU_CAN' HEU_presabs
# HI
awk '{print $1}' temp2 | grep -v "rep" | while read line; do if cat CAN_networks_hi_1_intersection.graphml | grep -w -m 1 -q $line; then
     echo "1"
else
     echo "0"
fi; 
done > HI_presabs
sed -i '1 i\HI_CAN' HI_presabs
# final dataframe
paste temp2 HUU_presabs HEU_presabs HI_presabs > cluster_groups.txt
```

Generate global cluster network with cluster colors and group margins

```R
set.seed(782469)
# install.packages("oaqc")
library(igraph)
library(ggforce)
library(concaveman)
library(ggraph)
library(oaqc)

dat <- read.table("full_CAN_all_networks_0.6_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
coords <- layout_nicely(can) 
# greedy community detection
clust <- cluster_fast_greedy(can)
# add in custom colors
# read in file that maps colors to global clusters -- colors added manually to this file!!
clustcol <- read.table("cluster_colors.txt", header=T)
# match node colors with attributes
V(can)$color <- clustcol$cluster_color[match(vertex_attr(can)$name, clustcol$repASV)]
V(can)$clust <- as.character(clustcol$cluster_simple[match(vertex_attr(can)$name, clustcol$repASV)])
# get a layout object
layout <- create_layout(can, layout="igraph", algorithm="nicely")

pdf("all_can_network.clust.pdf", width=10, height=10)
     ggraph(layout) +
     geom_mark_hull(aes(x=layout$x, y=layout$y, group=layout$clust, fill=layout$color, label=paste("Cluster ", layout$clust), filter = layout$clust != "other"), concavity=6, expand=unit(2, "mm"), alpha=0.25)+
     geom_edge_link0(aes(col = "grey"), width = 0.5) +
     geom_node_point(aes(color=as.factor(clust)), shape = 19, size = 7) +
     geom_node_point(aes(color="black"), shape = 21, size=7) +
     scale_color_manual(limits=as.factor(layout$clust), values=layout$color) +
     scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
     theme_graph()+
     theme(legend.position = "none")
dev.off()
```

Next, detect community modularity and clusters from our HIV group specific CAN networks

```bash
# convert graphml to edge list
grep "<edge " CAN_net_huu_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_huu_1_intersection.txt
grep "<edge " CAN_net_heu_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_heu_1_intersection.txt
grep "<edge " CAN_net_hi_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_hi_1_intersection.txt
```

Load into R and generate modularity graphs. What we want is a graph with the HIV group specific clusters indicated as the hull markers BUT the node color should indicate the original global cluster identity 

```R
set.seed(782469)
library(igraph)
library(ggforce)
library(concaveman)
library(ggraph)
library(oaqc)

# HUU CAN network
dat <- read.table("CAN_networks_huu_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.8022758

# get dataframe of representative ASVs with their cluster membership
clusts <- list()
for(i in 1:length(clust[])){
     col1 <- eval(parse(text=paste0("clust[]$`",i,"`", sep="")))
     reps <- rep(i, length(col1))
     clusts[[i]] <- cbind(col1, reps)
     }
cdat <- as.data.frame(do.call(rbind, clusts))
colnames(cdat) <- c("repASV", "cluster")

# add group specific clusters to network object
V(can)$clust_grp <- as.character(cdat$cluster[match(vertex_attr(can)$name, cdat$repASV)])
# read in global cluster color mapping file
clustcol <- read.table("cluster_colors.txt", header=T)
# add global cluster node colors to network object
V(can)$color <- clustcol$cluster_color[match(vertex_attr(can)$name, clustcol$repASV)]
# change nans to white
# replace NA with white
V(can)$color[is.na(V(can)$color)] <- "#FFFFFF"
# get a layout object
layout <- create_layout(can, layout="igraph", algorithm="nicely")

pdf("HUU_can_network.clust.pdf", width=10, height=10)
     ggraph(layout) +
     geom_mark_hull(aes(x=layout$x, y=layout$y, group=layout$clust_grp), concavity=10, expand=unit(2, "mm"), alpha=0.25)+
     geom_edge_link0(aes(col = "grey"), width = 0.5) +
     geom_node_point(aes(color=as.factor(clust_grp)), shape = 19, size = 7) +
     scale_color_manual(limits=as.factor(layout$clust_grp), values=layout$color) +
     scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
     geom_node_point(aes(color="black"), shape = 21, size=7) +
     theme_graph()+
     theme(legend.position = "none")
dev.off()
# write clusters to file
write.table(cdat, file="HUU_can_clusters.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

########################################
# HI CAN network
dat <- read.table("CAN_networks_hi_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.798711

# get dataframe of representative ASVs with their cluster membership
clusts <- list()
for(i in 1:length(clust[])){
     col1 <- eval(parse(text=paste0("clust[]$`",i,"`", sep="")))
     reps <- rep(i, length(col1))
     clusts[[i]] <- cbind(col1, reps)
     }
cdat <- as.data.frame(do.call(rbind, clusts))
colnames(cdat) <- c("repASV", "cluster")

# add group specific clusters to network object
V(can)$clust_grp <- as.character(cdat$cluster[match(vertex_attr(can)$name, cdat$repASV)])
# read in global cluster color mapping file
clustcol <- read.table("cluster_colors.txt", header=T)
# add global cluster node colors to network object
V(can)$color <- clustcol$cluster_color[match(vertex_attr(can)$name, clustcol$repASV)]
# change nans to white
# replace NA with white
V(can)$color[is.na(V(can)$color)] <- "#FFFFFF"
# get a layout object
layout <- create_layout(can, layout="igraph", algorithm="nicely")

pdf("HI_can_network.clust.pdf", width=10, height=10)
     ggraph(layout) +
     geom_mark_hull(aes(x=layout$x, y=layout$y, group=layout$clust_grp), concavity=10, expand=unit(2, "mm"), alpha=0.25)+
     geom_edge_link0(aes(col = "grey"), width = 0.5) +
     geom_node_point(aes(color=as.factor(clust_grp)), shape = 19, size = 7) +
     scale_color_manual(limits=as.factor(layout$clust_grp), values=layout$color) +
     scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
     geom_node_point(aes(color="black"), shape = 21, size=7) +
     theme_graph()+
     theme(legend.position = "none")
dev.off()
# write clusters to file
write.table(cdat, file="HI_can_clusters.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

########################################
# HEU CAN network
dat <- read.table("CAN_networks_heu_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.8312948

# get dataframe of representative ASVs with their cluster membership
clusts <- list()
for(i in 1:length(clust[])){
     col1 <- eval(parse(text=paste0("clust[]$`",i,"`", sep="")))
     reps <- rep(i, length(col1))
     clusts[[i]] <- cbind(col1, reps)
     }
cdat <- as.data.frame(do.call(rbind, clusts))
colnames(cdat) <- c("repASV", "cluster")

# add group specific clusters to network object
V(can)$clust_grp <- as.character(cdat$cluster[match(vertex_attr(can)$name, cdat$repASV)])
# read in global cluster color mapping file
clustcol <- read.table("cluster_colors.txt", header=T)
# add global cluster node colors to network object
V(can)$color <- clustcol$cluster_color[match(vertex_attr(can)$name, clustcol$repASV)]
# change nans to white
# replace NA with white
V(can)$color[is.na(V(can)$color)] <- "#FFFFFF"
# get a layout object
layout <- create_layout(can, layout="igraph", algorithm="nicely")

pdf("HEU_can_network.clust.pdf", width=10, height=10)
     ggraph(layout) +
     geom_mark_hull(aes(x=layout$x, y=layout$y, group=layout$clust_grp), concavity=10, expand=unit(2, "mm"), alpha=0.25)+
     geom_edge_link0(aes(col = "grey"), width = 0.5) +
     geom_node_point(aes(color=as.factor(clust_grp)), shape = 19, size = 7) +
     scale_color_manual(limits=as.factor(layout$clust_grp), values=layout$color) +
     scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
     geom_node_point(aes(color="black"), shape = 21, size=7) +
     theme_graph()+
     theme(legend.position = "none")
dev.off()
# write clusters to file
write.table(cdat, file="HEU_can_clusters.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
```

Healthy teeth in each HIV group

```bash
# convert graphml to edge list
grep "<edge " CAN_net_hi_H_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_hi_H_1_intersection.txt
grep "<edge " CAN_net_huu_H_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_huu_H_1_intersection.txt
grep "<edge " CAN_net_heu_H_1_intersection.graphml | sed 's/.*<edge source="//' | sed 's/" target="/\t/' | sed 's/">//' > CAN_networks_heu_H_1_intersection.txt
```

Load into R and generate modularity -- are there differences between groups?

```R
set.seed(782469)
library(igraph)
library(ggforce)
library(concaveman)
library(ggraph)
library(oaqc)

# HUU CAN network
dat <- read.table("CAN_networks_huu_H_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.206526

########################################
# HI CAN network
dat <- read.table("CAN_networks_hi_H_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.2012989

########################################
# HEU CAN network
dat <- read.table("CAN_networks_heu_H_1_intersection.txt", header=F)
can <- graph.data.frame(dat, directed=F)
# greedy community detection
clust <- cluster_fast_greedy(can)
modularity(clust)
# [1] 0.2160215
```
