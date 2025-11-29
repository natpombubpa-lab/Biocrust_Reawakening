#Title: Network inference using SpiecEasi package for uploading into HPC
#Authur: Teeratat Kaewjon; Pombubpa Lab
#Date: 6th October 2024

#load("./spiec_t500_n100_min91_LAC_c3_gl.rda")

library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)

library(ape)
library(vegan)
library(plyr)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)
library(tidyverse)
library(multcompView)
library(car)
library(igraph)
library(SpiecEasi)

#############################################################
####################### INPUT DATA  #########################
#############################################################
#Loading metabolite data
met.matrix <- read.table("./Metabolite_matrix.txt",header=T,sep="\t",row.names=1, check.names = FALSE)
metmat <- as(as.matrix(met.matrix), "matrix")
MET = otu_table(metmat, taxa_are_rows = TRUE)

taxmet <- read.table("./network.metabolite_Information.txt", header=T,sep="\t",row.names=1)
taxmet <- as(as.matrix(taxmet),"matrix")
TAX.met = tax_table(taxmet)

meta.met = read.table("./network.metabolite_metada.txt",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
sampleData.met <- sample_data(meta.met)

physeq.met = phyloseq(MET,TAX.met,sampleData.met)
physeq.met

met.log2 <- log2(otu_table(physeq.met))
met.log2

metmat.log2 <- as(as.matrix(met.log2), "matrix")
MET.log2 = otu_table(metmat.log2, taxa_are_rows = TRUE)
physeq.met = phyloseq(MET.log2,TAX.met,sampleData.met)
otu_table(physeq.met)

physeq.met <- subset_taxa(physeq.met, Compound != "Glutamine")
physeq.met <- subset_taxa(physeq.met, Compound != "Hypoxanthine")
physeq.met <- subset_taxa(physeq.met, Compound != "gamma-Aminobutyric acid")
physeq.met <- subset_taxa(physeq.met, Compound != "Glutamic acid")
physeq.met <- subset_taxa(physeq.met, Compound != "Methionine sulfoxide")
physeq.met <- subset_taxa(physeq.met, Compound != "Alanine")
physeq.met <- subset_taxa(physeq.met, Compound != "N,N-dimethylglycine")
physeq.met <- subset_taxa(physeq.met, Compound != "Hexoseamine")
physeq.met <- subset_taxa(physeq.met, Compound != "Aspartic acid")
physeq.met <- subset_taxa(physeq.met, Compound != "cyclic GMP")
physeq.met <- subset_taxa(physeq.met, Compound != "Thiamine")
physeq.met

#Loading ITS data

otus <- read.table("./CrustTemp2ITSasvUN3.otu_table_fix_header.txt",header=T,sep="\t",row.names=1)
otumat <- as(as.matrix(otus), "matrix")
OTU = otu_table(otumat, taxa_are_rows = TRUE)

taxmat <- read.table("./CrustTemp2ITSasvUN3.ASVs.taxonomy.fix.txt", header=T,sep="\t",row.names=1)
taxmat <- as(as.matrix(taxmat),"matrix")
TAX = tax_table(taxmat)

meta = read.table("./network.ITS_metada.txt",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
sampleData <- sample_data(meta)

physeq = phyloseq(OTU,TAX,sampleData)
physeq

physeq <- subset_taxa(physeq, Kingdom == "Fungi")
#physeq <- subset_taxa(physeq, Genus != "NA")
#physeq <- subset_taxa(physeq, Order == "Pleosporales")
physeq

physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune

#Loading 16S data
otus.bac <- read.table("./network.16S.asv.table_fix.txt", header = TRUE, sep = "\t", row.names = 1)
colnames(otus.bac) <- gsub("X", "", colnames(otus.bac))  # Remove "X" from all column names
otumat16s <- as.matrix(otus.bac)
OTU16s <- otu_table(otumat16s, taxa_are_rows = TRUE)

taxmat16s <- read.table("./taxonomy_split.txt", header=T,sep="\t",row.names=1)
taxmat16s <- as(as.matrix(taxmat16s),"matrix")
TAX16s <- tax_table(taxmat16s)

meta16s <- read.table("./network.ITS_metada.txt",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
sampleData16s <- sample_data(meta16s)

physeq16s <- phyloseq(OTU16s,TAX16s,sampleData16s)
physeq16s

physeq16s <- subset_taxa(physeq16s, Kingdom == "Bacteria")
physeq16s
physeq16s <- subset_taxa(physeq16s, Family != "Mitochondria")
physeq16s

#physeq16s <- subset_taxa(physeq16s, Species != "NA")
#physeq16s <- subset_taxa(physeq16s, Species != "uncultured bacterium")
#physeq16s <- subset_taxa(physeq16s, Species != "metagenome")

physeq16S.prune = prune_taxa(taxa_sums(physeq16s) > 1, physeq16s)
physeq16S.prune

###Subset samples by crust type (CLC)
physeqITS.LAC <- subset_samples(physeq.prune, Crusttype == "LAC")
physeq16S.LAC <- subset_samples(physeq16S.prune, Crusttype == "LAC")
Met.LAC <- subset_samples(physeq.met, Crusttype == "LAC")

# Subset by response Cluster 3
physeqITS.LAC.c3 <- subset_samples(physeqITS.LAC, Incubation == "96")

physeq16S.LAC.c3 <- subset_samples(physeq16S.LAC, Incubation == "96")

Met.LAC.c3 <- subset_samples(Met.LAC, Incubation == "96")

physeqITS.LAC.c3.top <- prune_taxa(names(sort(taxa_sums(physeqITS.LAC.c3), TRUE))[1:500], physeqITS.LAC.c3)
physeq16S.LAC.c3.top <- prune_taxa(names(sort(taxa_sums(physeq16S.LAC.c3), TRUE))[1:500], physeq16S.LAC.c3)

#############################################################
######### NETWORK CONSTRUCTION USING SPEIC-EASI  ############
#############################################################

parges <- list(rep.num = 20, seed = 10, ncores = 1)
spiec_LAC_c3_gl_4061 <- spiec.easi(list(physeqITS.LAC.c3.top, physeq16S.LAC.c3.top, Met.LAC.c3), method='glasso', nlambda=40,lambda.min.ratio= 6e-1, pulsar.params = parges)

#Save object
save(spiec_LAC_c3_gl_4061, file = "./20250227/spiec_LAC_c3_gl_4061.rda")

getStability(spiec_LAC_c3_gl_4061)
getOptInd(spiec_LAC_c3_gl_4061)
sum(getRefit(spiec_LAC_c3_gl_4061))/2
spiec_LAC_c3_gl_4061$select$stars$summary
spiec_LAC_c3_gl_4061$lambda

#############################################################
################## CREATE IGRAPH OBJECT  ####################
#############################################################

merge.Bac.Met.Fun <- merge_phyloseq(physeqITS.LAC.c3.top, physeq16S.LAC.c3.top,Met.LAC.c3)
merge.Bac.Met.Fun

ig.spiec_LAC_c3_gl_4061 <- adj2igraph(getRefit(spiec_LAC_c3_gl_4061),vertex.attr=list(name=taxa_names(merge.Bac.Met.Fun)),rmEmptyNodes = TRUE )
print(ig.spiec_LAC_c3_gl_4061)

#Save igraph object
save(ig.spiec_LAC_c3_gl_4061, file = "./20250227/ig.spiec_LAC_c3_gl_4091.rda")


#############################################################
############### PLOT NETWORK USING PHYLOSEQ  ################
#############################################################

plotnetwork =  plot_network(ig.spiec_LAC_c3_gl_4061, merge.Bac.Met.Fun,type = "taxa", shape = "Kingdom", color = "Kingdom")

pdf("./20250227/ig.spiec_LAC_c3_gl_4061.pdf", width = 10, height = 8)
plotnetwork
dev.off()


