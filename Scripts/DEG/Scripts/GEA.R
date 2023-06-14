if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("systemPipeR")
BiocManager::install("GO.db")

library(biomaRt)
library(systemPipeR)
library(GO.db)

setwd("C:/Users/dmjde/OneDrive/oud/Documenten/Dirk/School/Jaar 3/Minor/Project/DEG")

listMarts() # To choose the BioMart database


#For Plants (adjust for maize)
listMarts(host="https://plants.ensembl.org")
m <- useMart("plants_mart", host="https://plants.ensembl.org")
listDatasets(m)
m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="https://plants.ensembl.org")
listAttributes(m) # Choose the data types you want to download
go <- getBM(attributes=c("go_id", "ensembl_gene_id", "namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
go[1:4,]
dir.create("./GO")
write.table(go, "GO/GOannotationsBiomart_Athaliana.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
#write.table(go, "dif_SL_root_070719_k10_GO_IN.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


catdbAthaliana <- makeCATdb(myfile="GO/GOannotationsBiomart_Athaliana.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
save(catdbAthaliana, file="GO/catdbAthaliana.RData")
#Batch GO term enrichment analysis


#For Plants (adjust for maize)
load("GO/catdbAthaliana.RData")
#GeneID <- read.delim("GeneIDs.txt")
GeneID <- read.delim("dif_SL_root_070719_k10_GO_IN.txt")
# Cluster1 <- subset(GeneID, GeneID[2] == 1)
# Assuming the dataframe is named "GeneID"
colnames(GeneID)[1] <- "gene"
#subset_data <- GeneID[GeneID$k10cluster %in% c(1, 2, 9, 10), ]
Genes1 <- GeneID[GeneID$k10cluster == 1, ]
Cluster_1 <- Genes1[, 1]
Cluster_1 <- data.frame(Cluster_1)


Genes2 <- GeneID[GeneID$k10cluster == 2, ]
Cluster_2 <- Genes2[, 1]
Cluster_2 <- data.frame(Cluster_2)

Genes9 <- GeneID[GeneID$k10cluster == 9, ]
Cluster_9 <- Genes9[, 1]
Cluster_9 <- data.frame(Cluster_9)

Genes10 <- GeneID[GeneID$k10cluster == 10, ]
Cluster_10 <- Genes10[, 1]
Cluster_10 <- data.frame(Cluster_10)

# voor elk cluster apart de analyse doen, alleen de kolom gene ID selecteren.

#subset_data <- GeneID[GeneID$k10cluster == 1, ]
#BatchResult <- GOCluster_Report(catdb=catdbAthaliana, setlist=GeneID, id_type="gene", CLSZ=1, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResult1 <- GOCluster_Report(catdb=catdbAthaliana, setlist=Cluster_1, id_type="gene", CLSZ=1, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResult2 <- GOCluster_Report(catdb=catdbAthaliana, setlist=Cluster_2, id_type="gene", CLSZ=1, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResult9 <- GOCluster_Report(catdb=catdbAthaliana, setlist=Cluster_9, id_type="gene", CLSZ=1, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResult10 <- GOCluster_Report(catdb=catdbAthaliana, setlist=Cluster_10, id_type="gene", CLSZ=1, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)

write.table(BatchResult1, "GO/BatchResultsGeneIDs_1.xls", quote=FALSE, sep="\t", col.names = NA)
write.table(BatchResult2, "GO/BatchResultsGeneIDs_2.xls", quote=FALSE, sep="\t", col.names = NA)
write.table(BatchResult9, "GO/BatchResultsGeneIDs_9.xls", quote=FALSE, sep="\t", col.names = NA)
write.table(BatchResult10, "GO/BatchResultsGeneIDs_10.xls", quote=FALSE, sep="\t", col.names = NA)

#Plot batch GO term results for categories with a p-value of <=0.05

tussen1 <- BatchResult1
gos1 <- tussen1[tussen1$Padj <= 0.05, ]
gos1 <- na.omit(gos1)

pdf("GOslimbarplotMF_1.pdf", height=8, width=10); 
goBarplot(gos1, gocat="MF"); dev.off()
pdf("GOslimbarplotBP_1.pdf", height=8, width=10); 
goBarplot(gos1, gocat="BP"); dev.off()
pdf("GOslimbarplotCC_1.pdf", height=8, width=10); 
goBarplot(gos1, gocat="CC"); dev.off()


tussen2 <- BatchResult2
gos2 <- tussen2[tussen2$Padj <= 0.05, ]
gos2 <- na.omit(gos2)

pdf("GOslimbarplotMF_2.pdf", height=8, width=10); 
goBarplot(gos2, gocat="MF"); dev.off()
pdf("GOslimbarplotBP_2.pdf", height=8, width=10); 
goBarplot(gos2, gocat="BP"); dev.off()
pdf("GOslimbarplotCC_2.pdf", height=8, width=10); 
goBarplot(gos2, gocat="CC"); dev.off()


tussen9 <- BatchResult9
gos9 <- tussen9[tussen9$Padj <= 0.05, ]
gos9 <- na.omit(gos9)

pdf("GOslimbarplotMF_9.pdf", height=8, width=10); 
goBarplot(gos9, gocat="MF"); dev.off()
pdf("GOslimbarplotBP_9.pdf", height=8, width=10); 
goBarplot(gos9, gocat="BP"); dev.off()
pdf("GOslimbarplotCC_9.pdf", height=8, width=10); 
goBarplot(gos9, gocat="CC"); dev.off()

tussen10 <- BatchResult10
gos10 <- tussen10[tussen10$Padj <= 0.05, ]
gos10 <- na.omit(gos10)

pdf("GOslimbarplotMF_10.pdf", height=8, width=10); 
goBarplot(gos10, gocat="MF"); dev.off()
pdf("GOslimbarplotBP_10.pdf", height=8, width=10); 
goBarplot(gos10, gocat="BP"); dev.off()
pdf("GOslimbarplotCC_10.pdf", height=8, width=10); 
goBarplot(gos10, gocat="CC"); dev.off()
