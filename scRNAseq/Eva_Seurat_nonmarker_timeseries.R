library(Seurat)
library(ggplot2)
library(cowplot)

SIGAH12.data <- Read10X(data.dir = "/Users/christopherpenfold/Desktop/10X/Eva_10X/SIGAH12_91_Chm3_FC2/outs/filtered_feature_bc_matrix")
SIGAH12 <- CreateSeuratObject(counts = SIGAH12.data, project = "SIGAH12")
SIGAH12$stim <- "SIGAH12"
SIGAH12 <- subset(SIGAH12, subset = nFeature_RNA > 0)
SIGAH12 <- NormalizeData(SIGAH12, verbose = FALSE)
SIGAH12 <- FindVariableFeatures(SIGAH12, selection.method = "vst", nfeatures = 2000)

SIGAH5.data <- Read10X(data.dir = "/Users/christopherpenfold/Desktop/10X/Eva_10X/SIGAH5_91_Chm3_FC2/outs/filtered_feature_bc_matrix")
SIGAH5 <- CreateSeuratObject(counts = SIGAH5.data, project = "SIGAH5")
SIGAH5$stim <- "SIGAH5"
SIGAH5 <- subset(SIGAH5, subset = nFeature_RNA > 0)
SIGAH5 <- NormalizeData(SIGAH5, verbose = FALSE)
SIGAH5 <- FindVariableFeatures(SIGAH5, selection.method = "vst", nfeatures = 2000)

#SIGAH12.data <- Read10X(data.dir = "/Users/christopherpenfold/Desktop/Eva/10X/SIGAH12_91_Chm3_FCr/outs/filtered_feature_bc_matrix")
#SIGAC11 <- CreateSeuratObject(counts = SIGAC11.data, project = "SIGAC11")
#SIGAC11$stim <- "SIGAC11"
#SIGAC11$isnt <- "IVF"
#SIGAC11 <- NormalizeData(SIGAC11, verbose = FALSE)
##SIGAC11 <- subset(SIGAC11, subset = nFeature_RNA > 500)
#SIGAC11 <- FindVariableFeatures(SIGAC11, selection.method = "vst", nfeatures = 2000)


#SIGAC11.data <- Read10X(data.dir = "/Users/christopherpenfold/Desktop/Eva/10X/SIGAC11_91_Chm3_FCr/outs/filtered_feature_bc_matrix")
#SIGAC11$stim <- "SIGAC11"
#SIGAC11 <- CreateSeuratObject(counts = SIGAC11.data, project = "SIGAC11")
#SIGAC11 <- subset(SIGAC11, subset = nFeature_RNA > 500)
#SIGAC11 <- NormalizeData(SIGAC11, verbose = FALSE)
#SIGAC11 <- FindVariableFeatures(SIGAC11, selection.method = "vst", nfeatures = 2000)

#SIGAD11.data <- Read10X(data.dir = "/Users/christopherpenfold/Desktop/Eva/10X/SIGAD11_91_Chm3_FCr/outs/filtered_feature_bc_matrix")
#SIGAD11 <- CreateSeuratObject(counts = SIGAD11.data, project = "SIGAC11")
#SIGAD11$stim <- "SIGAD11"
#SIGAD11 <- NormalizeData(SIGAD11, verbose = FALSE)
#SIGAD11 <- FindVariableFeatures(SIGAD11, selection.method = "vst", nfeatures = 2000)
#SIGAD11 <- subset(SIGAD11, subset = nFeature_RNA > 500)

Markers <- read.delim("/Users/christopherpenfold/Desktop/Eva/Markers.txt")
markerlist <- c(paste(Markers$Marker_genes,".L",sep = ""),paste(Markers$Marker_genes,".S",sep = ""))

Markers <- read.delim("/Users/christopherpenfold/Desktop/Eva/GeneExpression_laevis.txt", header = 0)


#gene_List2 <- rownames(SIGAH12)
#gene_List1 <- rownames(SIGAH5)
#anchorlist <- intersect(intersect(gene_List1,gene_List2),markerlist)

xenopus.anchors <- FindIntegrationAnchors(object.list = list(SIGAH5, SIGAH12), dims = 1:20, anchor.features = 2000)
xenopus.combined <- IntegrateData(anchorset = xenopus.anchors, dims = 1:20)
DefaultAssay(xenopus.combined) <- "integrated"
xenopus.combined <- ScaleData(xenopus.combined, verbose = FALSE)
xenopus.combined <- RunPCA(xenopus.combined, npcs = 30, verbose = FALSE)
xenopus.combined <- RunUMAP(xenopus.combined, reduction = "pca", dims = 1:20)
xenopus.combined <- RunTSNE(xenopus.combined, reduction = "pca", dims = 1:20)
#xenopus.combined <- RunDiffusion(xenopus.combined, reduction = "pca", dims = 1:20)

#Clustering etc
xenopus.combined <- FindNeighbors(xenopus.combined, reduction = "pca", dims = 1:20)
xenopus.combined <- FindClusters(xenopus.combined, resolution = 0.5)

p1 <- DimPlot(xenopus.combined, reduction = "umap", group.by = "stim", dim.1 = 10, dim.2 = 10)
p2 <- DimPlot(xenopus.combined, reduction = "umap", label = TRUE, dim.1 = 10, dim.2 = 10)
plt <- plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","UMAP",".pdf",sep=""),plot = plt)
#save_plot(filename=paste("~/Desktop/Eva/","UMAP",".pdf",sep=""),plt, ncol = 1, nrow = 2, base_height = 4,
#          base_aspect_ratio = 1.1)
dev.off()  


#p1 <- DimPlot(xenopus.combined, reduction = "dm", group.by = "stim", dim.1 = 10, dim.2 = 10)
#p2 <- DimPlot(xenopus.combined, reduction = "dm", label = TRUE, dim.1 = 10, dim.2 = 10)
#plt <- plot_grid(p1, p2)
#ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DM",".pdf",sep=""),plot = plt)
#save_plot(filename=paste("~/Desktop/Eva/","UMAP",".pdf",sep=""),plt, ncol = 1, nrow = 2, base_height = 4,
#          base_aspect_ratio = 1.1)
#dev.off()  


#xenopus.anchors <- FindIntegrationAnchors(object.list = list(SIGAH12, SIGAH5), dims = 1:20, anchor.features = 2000)
#xenopus.combined <- IntegrateData(anchorset = xenopus.anchors, dims = 1:20)

#DefaultAssay(xenopus.combined) <- "integrated"

#Generate the reduced dimensionality coords
#xenopus.combined <- ScaleData(xenopus.combined, verbose = FALSE)
#xenopus.combined <- RunPCA(xenopus.combined, npcs = 30, verbose = FALSE)
#xenopus.combined <- RunUMAP(xenopus.combined, reduction = "pca", dims = 1:20)
#xenopus.combined <- RunTSNE(xenopus.combined, reduction = "pca", dims = 1:20)
#Clustering etc
#xenopus.combined <- FindNeighbors(xenopus.combined, reduction = "pca", dims = 1:20)
#xenopus.combined <- FindClusters(xenopus.combined, resolution = 0.5)

#p1 <- DimPlot(xenopus.combined, reduction = "umap", group.by = "stim", dim.1 = 10, dim.2 = 10)
#p2 <- DimPlot(xenopus.combined, reduction = "umap", label = TRUE, dim.1 = 10, dim.2 = 10)
#plt <- plot_grid(p1, p2)
#ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","UMAP",".pdf",sep=""),plot = plt)
#save_plot(filename=paste("~/Desktop/Eva/","UMAP",".pdf",sep=""),plt, ncol = 1, nrow = 2, base_height = 4,
#          base_aspect_ratio = 1.1)
#dev.off()          

write.csv(as.data.frame(Idents(object = xenopus.combined)), file="~/Desktop/Eva/OnClusterPlots/EmbeddingsCl.csv")
write.csv(as.data.frame(Embeddings(object = xenopus.combined[["umap"]])), file="~/Desktop/Eva/OnClusterPlots/Embeddings.csv")
#write.csv(as.data.frame(Embeddings(object = xenopus.combined[["dm"]])), file="~/Desktop/Eva/OnClusterPlots/EmbeddingsDM.csv")
write.csv(as.data.frame(Embeddings(object = xenopus.combined[["tsne"]])), file="~/Desktop/Eva/OnClusterPlots/EmbeddingsTSNE.csv")
write.csv(as.data.frame(xenopus.combined[[]]), file="~/Desktop/Eva/OnClusterPlots/EmbeddingsKey.csv")

p1 <- DimPlot(xenopus.combined, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(xenopus.combined, reduction = "tsne", label = TRUE)
plt <- plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","TSNE_Split",".pdf",sep=""),plot = plt)
dev.off()

DimPlot(xenopus.combined, reduction = "umap", split.by = "stim")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","UMAP_Split",".pdf",sep=""))
dev.off() 


DimPlot(xenopus.combined, reduction = "dm", split.by = "stim")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DM_Split",".pdf",sep=""))
dev.off() 

DimPlot(xenopus.combined, reduction = "tsne", split.by = "stim")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","TSNE_Split",".pdf",sep=""))
dev.off() 


#
early = c("animal cap","dorsal marginal zone","neuroectoderm","non-neural ectoderm","marginal zone","endoderm","foregut","hindgut","organizer","notochord","involuted ventral mesoderm","involuted dorsal mesoderm","tail bud","cement gland primordium","goblet cell","ionocyte","ciliated epidermal progenitor","neural plate","neural plate border")

OnMem <- read.delim("/Users/christopherpenfold/Desktop/Eva/Motifs/Genelists_from_Eva_April2018/EndoNT_Filter_DonorsvsIVFdownFDR_DonsRPKM1_IVFvsNT_FDRFC3up__OnMem3FC.txt")
OffMem <- read.delim("/Users/christopherpenfold/Desktop/Eva/Motifs/Genelists_from_Eva_April2018/EndoNT_Filter_DonorsvsIVFupFDR_IVFvsNT_FDRFC3down__OffMem3FC.txt")


mat1 <- as.matrix(GetAssayData(SIGAH12, slot = "counts"))
N_1 <- matrix(0, nrow = dim(mat1)[2], ncol = length(early)+2)
N_2 <- matrix(0, nrow = dim(mat1)[2], ncol = length(early)+2)
N_3 <- matrix(0, nrow = dim(mat1)[2], ncol = length(early)+2)
for (j in 1:length(early)){
  mg <- Markers$V2[grep(as.character(early[j]),as.character(Markers$V4))]
    for (i in 1:dim(mat1)[2]){
      N_1[i,j] <- length(intersect(rownames(mat1)[which(mat1[,i]>5)],mg))
      N_2[i,j] <- length((rownames(mat1)[which(mat1[,i]>5)]))      
      ps <- phyper(N_1[i,j]-1, N_2[i,j], dim(mat1)[1]-N_2[i,j], length(mg), lower.tail= FALSE)
      N_3[i,j] <- ps
    }
}

for (i in 1:dim(mat1)[2]){
N_1[i,j+1] <- length(intersect(rownames(mat1)[which(mat1[,i]>5)],OnMem[,1]))
N_2[i,j+1] <- length((rownames(mat1)[which(mat1[,i]>5)]))      
ps <- phyper(N_1[i,j+1]-1, N_2[i,j+1], dim(mat1)[1]-N_2[i,j+1], dim(OnMem)[1], lower.tail= FALSE)
N_3[i,j+1] <- ps
}
for (i in 1:dim(mat1)[2]){
  N_1[i,j+2] <- length(intersect(rownames(mat1)[which(mat1[,i]>5)],OffMem[,1]))
  N_2[i,j+2] <- length((rownames(mat1)[which(mat1[,i]>5)]))      
  ps <- phyper(N_1[i,j+2]-1, N_2[i,j+2], dim(mat1)[1]-N_2[i,j+2], dim(OffMem)[1], lower.tail= FALSE)
  N_3[i,j+2] <- ps
}



mat2 <- as.matrix(GetAssayData(SIGAH5, slot = "counts"))
M_1 <- matrix(0, nrow = dim(mat2)[2], ncol = length(early)+2)
M_2 <- matrix(0, nrow = dim(mat2)[2], ncol = length(early)+2)
M_3 <- matrix(0, nrow = dim(mat2)[2], ncol = length(early)+2)
for (j in 1:length(early)){
  mg <- Markers$V2[grep(as.character(early[j]),as.character(Markers$V4))]
  for (i in 1:dim(mat2)[2]){
    M_1[i,j] <- length(intersect(rownames(mat2)[which(mat2[,i]>5)],mg))
    M_2[i,j] <- length((rownames(mat2)[which(mat2[,i]>5)]))      
    #M_3[i,j] <- length((mg))
    ps <- phyper(M_1[i,j]-1, M_2[i,j], dim(mat2)[1]-M_2[i,j], length(mg), lower.tail= FALSE)
    M_3[i,j] <- ps
    
    
  }
}

for (i in 1:dim(mat1)[2]){
  M_1[i,j+1] <- length(intersect(rownames(mat2)[which(mat2[,i]>5)],OnMem[,1]))
  M_2[i,j+1] <- length((rownames(mat2)[which(mat2[,i]>5)]))      
  ps <- phyper(M_1[i,j+1]-1, M_2[i,j+1], dim(mat2)[1]-M_2[i,j+1], dim(OnMem)[1], lower.tail= FALSE)
  M_3[i,j+1] <- ps
}
for (i in 1:dim(mat1)[2]){
  M_1[i,j+2] <- length(intersect(rownames(mat2)[which(mat2[,i]>5)],OffMem[,1]))
  M_2[i,j+2] <- length((rownames(mat2)[which(mat2[,i]>5)]))      
  ps <- phyper(M_1[i,j+2]-1, M_2[i,j+2], dim(mat2)[1]-M_2[i,j+2], dim(OffMem)[1], lower.tail= FALSE)
  M_3[i,j+2] <- ps
}

early = c("animal cap","dorsal marginal zone","neuroectoderm","non-neural ectoderm","marginal zone","endoderm","foregut","hindgut","organizer","notochord","involuted ventral mesoderm","involuted dorsal mesoderm","tail bud","cement gland primordium","goblet cell","ionocyte","ciliated epidermal progenitor","neural plate","neural plate border")


early <- c(early,"On","Off")

for (j in 1:length(early)){

  pvalss <- c(N_3[,j],M_3[,j])
  padj <- pvalss #p.adjust(pvalss, method = "bonferroni")
  
  vec2add <- 1-padj #c(N_3[,j],M_3[,j]) #as.data.frame(t(c(N_3[,j],M_3[,j])))
names(vec2add) <- colnames(x = xenopus.combined)
#cluster_letters <- LETTERS[Idents(object = xenopus.combined)]
#names(cluster_letters) <- colnames(x = xenopus.combined)
xenopus.combined <- AddMetaData(
  object = xenopus.combined,
  metadata = vec2add,
  col.name = 'stuff'
)
FeaturePlot(xenopus.combined, features = 'stuff', split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/",early[j],j,".pdf",sep=""))
}

for (j in 1:length(early)){
  vec2add <- as.data.frame(c(N_3[,j],M_3[,j]))
  rownames(vec2add) <- colnames(xenopus.combined)
AddMetaData(xenopus.combined, vec2add, col.name = early[j])
FeaturePlot(xenopus.combined, features = early[j], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
}  

for (j in 1:length(early)){
FeaturePlot(xenopus.combined, features = early[1], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/",early[i], "_Overrep_",i,".pdf",sep=""))
}



#Look at cluster 7 markers vs others and plot the top X
DefaultAssay(xenopus.combined) <- "RNA"

cl0.markers <- FindMarkers(xenopus.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
cl1.markers <- FindMarkers(xenopus.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
cl2.markers <- FindMarkers(xenopus.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
cl3.markers <- FindMarkers(xenopus.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
cl4.markers <- FindMarkers(xenopus.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
cl5.markers <- FindMarkers(xenopus.combined, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
cl6.markers <- FindMarkers(xenopus.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
cl7.markers <- FindMarkers(xenopus.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
cl8.markers <- FindMarkers(xenopus.combined, ident.1 = 8, grouping.var = "stim", verbose = FALSE)
cl9.markers <- FindMarkers(xenopus.combined, ident.1 = 9, grouping.var = "stim", verbose = FALSE)

write.csv(as.data.frame(cl0.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_0.csv",row.names=FALSE)
write.csv(as.data.frame(cl1.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_1.csv",row.names=FALSE)
write.csv(as.data.frame(cl2.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_2.csv",row.names=FALSE)
write.csv(as.data.frame(cl3.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_3.csv",row.names=FALSE)
write.csv(as.data.frame(cl4.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_4.csv",row.names=FALSE)
write.csv(as.data.frame(cl5.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_5.csv",row.names=FALSE)
write.csv(as.data.frame(cl6.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_6.csv",row.names=FALSE)
write.csv(as.data.frame(cl7.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_7.csv",row.names=FALSE)
write.csv(as.data.frame(cl8.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_8.csv",row.names=FALSE)
write.csv(as.data.frame(cl9.markers), file = "~/Desktop/Eva/OnClusterPlots/ClDEmarker_9.csv",row.names=FALSE)

cl0.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
cl1.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
cl2.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
cl3.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
cl4.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
cl5.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
cl6.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
cl7.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
cl8.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 8, grouping.var = "stim", verbose = FALSE)
cl9.cmarkers <- FindConservedMarkers(xenopus.combined, ident.1 = 9, grouping.var = "stim", verbose = FALSE)
#cl10.markers <- FindConservedMarkers(xenopus.combined, ident.1 = 10, grouping.var = "stim", verbose = FALSE)
#cl11.markers <- FindConservedMarkers(xenopus.combined, ident.1 = 11, grouping.var = "stim", verbose = FALSE)
write.csv(as.data.frame(cl0.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_0.csv",row.names=FALSE)
write.csv(as.data.frame(cl1.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_1.csv",row.names=FALSE)
write.csv(as.data.frame(cl2.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_2.csv",row.names=FALSE)
write.csv(as.data.frame(cl3.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_3.csv",row.names=FALSE)
write.csv(as.data.frame(cl4.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_4.csv",row.names=FALSE)
write.csv(as.data.frame(cl5.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_5.csv",row.names=FALSE)
write.csv(as.data.frame(cl6.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_6.csv",row.names=FALSE)
write.csv(as.data.frame(cl7.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_7.csv",row.names=FALSE)
write.csv(as.data.frame(cl8.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_8.csv",row.names=FALSE)
write.csv(as.data.frame(cl9.cmarkers), file = "~/Desktop/Eva/OnClusterPlots/ClConmarker_9.csv",row.names=FALSE)




FeaturePlot(xenopus.combined, features = rownames(cl0.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl0",".pdf",sep=""))

FeaturePlot(xenopus.combined, features = rownames(cl1.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl1",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl2.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl2",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl3.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl3",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl4.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl4",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl5.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl5",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl6.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl6",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl7.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl7",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl8.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl8",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl9.markers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl9",".pdf",sep=""))
#FeaturePlot(xenopus.combined, features = rownames(cl10.markers)[1:4], min.cutoff = "q9")
#ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl10",".pdf",sep=""))


FeaturePlot(xenopus.combined, features = rownames(cl0.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl0",".pdf",sep=""))

FeaturePlot(xenopus.combined, features = rownames(cl1.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl1",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl2.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl2",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl3.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl3",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl4.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl4",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl5.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl5",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl6.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl6",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl7.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl7",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl8.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl8",".pdf",sep=""))
FeaturePlot(xenopus.combined, features = rownames(cl9.cmarkers)[1:4], min.cutoff = "q9")
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","cCl9",".pdf",sep=""))



for (i in 1:200){
  FeaturePlot(xenopus.combined, features = rownames(cl1.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl1__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl1.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl1__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl2.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl2__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl3.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl3__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl4.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl4__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl5.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl5__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl6.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl6__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl7.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl7__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl8.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl8__",i,".pdf",sep=""))
  FeaturePlot(xenopus.combined, features = rownames(cl9.markers)[i], min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl9__",i,".pdf",sep=""))
  #FeaturePlot(xenopus.combined, features = rownames(cl10.markers)[i], min.cutoff = "q9")
  ##ggsave(filename=paste("~/Desktop/Eva/","Cl10__",i,".pdf",sep=""))}
}



markers.to.plot = unique(c(rownames(cl0.markers)[1:5],
                    rownames(cl1.markers)[1:5],
                    rownames(cl2.markers)[1:5],
                    rownames(cl3.markers)[1:5],
                    rownames(cl4.markers)[1:5],
                    rownames(cl5.markers)[1:5],
                    rownames(cl6.markers)[1:5],
                    rownames(cl7.markers)[1:5],
                    rownames(cl8.markers)[1:5],
                    rownames(cl9.markers)[1:5]
                    ))

DotPlot(xenopus.combined, features = rev(markers.to.plot), cols = c("blue", "red","green","yellow"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","C0-l0",".pdf",sep=""))



markers.to.plot = unique(c(rownames(cl0.markers)[6:10],
                           rownames(cl1.markers)[6:10],
                           rownames(cl2.markers)[6:10],
                           rownames(cl3.markers)[6:10],
                           rownames(cl4.markers)[6:10],
                           rownames(cl5.markers)[6:10],
                           rownames(cl6.markers)[6:10],
                           rownames(cl7.markers)[6:10],
                           rownames(cl8.markers)[6:10],
                           rownames(cl9.markers)[6:10]
))

DotPlot(xenopus.combined, features = rev(markers.to.plot), cols = c("blue", "red","green","yellow"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","C0-l0_2",".pdf",sep=""))


markers.to.plot = unique(c(rownames(cl0.markers)[11:15],
                           rownames(cl1.markers)[11:15],
                           rownames(cl2.markers)[11:15],
                           rownames(cl3.markers)[11:15],
                           rownames(cl4.markers)[11:15],
                           rownames(cl5.markers)[11:15],
                           rownames(cl6.markers)[11:15],
                           rownames(cl7.markers)[11:15],
                           rownames(cl8.markers)[11:15],
                           rownames(cl9.markers)[11:15]
))

DotPlot(xenopus.combined, features = rev(markers.to.plot), cols = c("blue", "red","green","yellow"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","C0-l0_3",".pdf",sep=""))


markers.to.plot = unique(c(rownames(cl0.markers)[16:20],
                           rownames(cl1.markers)[16:20],
                           rownames(cl2.markers)[16:20],
                           rownames(cl3.markers)[16:20],
                           rownames(cl4.markers)[16:20],
                           rownames(cl5.markers)[16:20],
                           rownames(cl6.markers)[16:20],
                           rownames(cl7.markers)[16:20],
                           rownames(cl8.markers)[16:20],
                           rownames(cl9.markers)[16:20]
))

DotPlot(xenopus.combined, features = rev(markers.to.plot), cols = c("blue", "red","green","yellow"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","C0-l0_4",".pdf",sep=""))


#immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
#                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
#                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")
#
#DimPlot(immune.combined, label = TRUE)


#
#Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", 
#                                                                      "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
#markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
#                                                                      "CD4 Naive T", "CD4 Memory T"))
#                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
#                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
#DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
#        split.by = "stim") + RotatedAxis()


cl1.cells <- subset(xenopus.combined, idents = "1")
Idents(cl1.cells) <- "stim"
avg.cl1.cells <- log1p(AverageExpression(cl1.cells, verbose = FALSE)$RNA)
inds1 = rownames(avg.cl1.cells)[which(avg.cl1.cells$SIGAH12-avg.cl1.cells$SIGAH5> log(2) | avg.cl1.cells$SIGAH5-avg.cl1.cells$SIGAH12 > log(2))]
avg.cl1.cells$gene <- rownames(avg.cl1.cells)

cl2.cells <- subset(xenopus.combined, idents = "2")
Idents(cl2.cells) <- "stim"
avg.cl2.cells <- log1p(AverageExpression(cl2.cells, verbose = FALSE)$RNA)
inds2 = rownames(avg.cl2.cells)[which(avg.cl2.cells$SIGAH12-avg.cl2.cells$SIGAH5> log(2) | avg.cl2.cells$SIGAH5-avg.cl2.cells$SIGAH12 > log(2))]
avg.cl2.cells$gene <- rownames(avg.cl2.cells)

#genes.to.label = c("Xelaev18000048m.g", "pga4.L", "hes6.1.S", "cdx4.L", "hoxd1.L", "cbfb.S")
p1 <- ggplot(avg.cl1.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl1 Cells")
p2 <- ggplot(avg.cl2.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl2 Cells")
p1 <- LabelPoints(plot = p1, points = inds1, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = inds2, repel = TRUE)
plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_Cl1_Cl2",".pdf",sep=""))


cl3.cells <- subset(xenopus.combined, idents = "3")
Idents(cl3.cells) <- "stim"
avg.cl3.cells <- log1p(AverageExpression(cl3.cells, verbose = FALSE)$RNA)
inds3 = rownames(avg.cl3.cells)[which(avg.cl3.cells$SIGAH12-avg.cl3.cells$SIGAH5> log(2) | avg.cl3.cells$SIGAH5-avg.cl3.cells$SIGAH12 > log(2))]
avg.cl3.cells$gene <- rownames(avg.cl3.cells)

cl4.cells <- subset(xenopus.combined, idents = "4")
Idents(cl4.cells) <- "stim"
avg.cl4.cells <- log1p(AverageExpression(cl4.cells, verbose = FALSE)$RNA)
inds4 = rownames(avg.cl4.cells)[which(avg.cl4.cells$SIGAH12-avg.cl4.cells$SIGAH5> log(2) | avg.cl4.cells$SIGAH5-avg.cl4.cells$SIGAH12 > log(2))]
avg.cl4.cells$gene <- rownames(avg.cl4.cells)

#genes.to.label = c("Xelaev18000048m.g", "pga4.L", "hes6.1.S", "cdx4.L", "hoxd1.L", "cbfb.S")
p1 <- ggplot(avg.cl3.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl3 Cells")
p2 <- ggplot(avg.cl4.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl4 Cells")
p1 <- LabelPoints(plot = p1, points = inds3, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = inds4, repel = TRUE)
plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_Cl3_Cl4",".pdf",sep=""))



cl5.cells <- subset(xenopus.combined, idents = "5")
Idents(cl5.cells) <- "stim"
avg.cl5.cells <- log1p(AverageExpression(cl5.cells, verbose = FALSE)$RNA)
inds5 = rownames(avg.cl5.cells)[which(avg.cl5.cells$SIGAH12-avg.cl5.cells$SIGAH5> log(2) | avg.cl5.cells$SIGAH5-avg.cl5.cells$SIGAH12 > log(2))]
avg.cl5.cells$gene <- rownames(avg.cl5.cells)

cl6.cells <- subset(xenopus.combined, idents = "6")
Idents(cl6.cells) <- "stim"
avg.cl6.cells <- log1p(AverageExpression(cl6.cells, verbose = FALSE)$RNA)
inds6 = rownames(avg.cl6.cells)[which(avg.cl6.cells$SIGAH12-avg.cl6.cells$SIGAH5> log(2) | avg.cl6.cells$SIGAH5-avg.cl6.cells$SIGAH12 > log(2))]
avg.cl6.cells$gene <- rownames(avg.cl6.cells)

#genes.to.label = c("Xelaev18000048m.g", "pga4.L", "hes6.1.S", "cdx4.L", "hoxd1.L", "cbfb.S")
p1 <- ggplot(avg.cl5.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl5 Cells")
p2 <- ggplot(avg.cl6.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl6 Cells")
p1 <- LabelPoints(plot = p1, points = inds5, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = inds6, repel = TRUE)
plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_Cl5_Cl6",".pdf",sep=""))



cl7.cells <- subset(xenopus.combined, idents = "7")
Idents(cl7.cells) <- "stim"
avg.cl7.cells <- log1p(AverageExpression(cl7.cells, verbose = FALSE)$RNA)
inds7 = rownames(avg.cl7.cells)[which(avg.cl7.cells$SIGAH12-avg.cl7.cells$SIGAH5> log(2) | avg.cl7.cells$SIGAH5-avg.cl7.cells$SIGAH12 > log(2))]
avg.cl7.cells$gene <- rownames(avg.cl7.cells)

cl8.cells <- subset(xenopus.combined, idents = "8")
Idents(cl8.cells) <- "stim"
avg.cl8.cells <- log1p(AverageExpression(cl8.cells, verbose = FALSE)$RNA)
inds8 = rownames(avg.cl8.cells)[which(avg.cl8.cells$SIGAH12-avg.cl8.cells$SIGAH5> log(2) | avg.cl8.cells$SIGAH5-avg.cl8.cells$SIGAH12 > log(2))]
avg.cl8.cells$gene <- rownames(avg.cl8.cells)

#genes.to.label = c("Xelaev18000048m.g", "pga4.L", "hes6.1.S", "cdx4.L", "hoxd1.L", "cbfb.S")
p1 <- ggplot(avg.cl7.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl5 Cells")
p2 <- ggplot(avg.cl8.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl6 Cells")
p1 <- LabelPoints(plot = p1, points = inds7, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = inds8, repel = TRUE)
plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_Cl7_Cl8",".pdf",sep=""))



cl9.cells <- subset(xenopus.combined, idents = "9")
Idents(cl9.cells) <- "stim"
avg.cl9.cells <- log1p(AverageExpression(cl9.cells, verbose = FALSE)$RNA)
inds9 = rownames(avg.cl9.cells)[which(avg.cl9.cells$SIGAH12-avg.cl9.cells$SIGAH5> log(2) | avg.cl9.cells$SIGAH5-avg.cl9.cells$SIGAH12 > log(2))]
avg.cl9.cells$gene <- rownames(avg.cl9.cells)

#cl10.cells <- subset(xenopus.combined, idents = "10")
#Idents(cl10.cells) <- "stim"
#avg.cl10.cells <- log1p(AverageExpression(cl10.cells, verbose = FALSE)$RNA)
#inds10 = rownames(avg.cl10.cells)[which(avg.cl10.cells$SIGAH12-avg.cl10.cells$SIGAH5> log(2) | avg.cl10.cells$SIGAH5-avg.cl8.cells$SIGAH12 > log(2))]
#avg.cl10.cells$gene <- rownames(avg.cl10.cells)

#genes.to.label = c("Xelaev18000048m.g", "pga4.L", "hes6.1.S", "cdx4.L", "hoxd1.L", "cbfb.S")
#p1 <- ggplot(avg.cl9.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl5 Cells")
#p2 <- ggplot(avg.cl10.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl6 Cells")
#p1 <- LabelPoints(plot = p1, points = inds9, repel = TRUE)
#p2 <- LabelPoints(plot = p2, points = inds10, repel = TRUE)
#plot_grid(p1, p2)
#ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_Cl9_Cl10",".pdf",sep=""))

cl0.cells <- subset(xenopus.combined, idents = "9")
Idents(cl0.cells) <- "stim"
avg.cl0.cells <- log1p(AverageExpression(cl0.cells, verbose = FALSE)$RNA)
inds0 = rownames(avg.cl0.cells)[which(avg.cl0.cells$SIGAH12-avg.cl0.cells$SIGAH5> log(2) | avg.cl0.cells$SIGAH5-avg.cl0.cells$SIGAH12 > log(2))]
avg.cl0.cells$gene <- rownames(avg.cl0.cells)

p1 <- ggplot(avg.cl9.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl9 Cells")
p2 <- ggplot(avg.cl0.cells, aes(SIGAH12, SIGAH5)) + geom_point() + ggtitle("Cl0 Cells")
p1 <- LabelPoints(plot = p1, points = inds9, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = inds0, repel = TRUE)
#plot_grid(p1, p2)
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","DE_C9_Cl0",".pdf",sep=""))


DElist = c(inds0,inds1,inds2,inds3,inds4,inds5,inds6,inds7,inds8,inds9)
labels =  c(rep('0',length(inds0)),rep('1',length(inds1)),
            rep('2',length(inds2)),rep('3',length(inds3)),
            rep('4',length(inds4)),rep('5',length(inds5)),
            rep('6',length(inds6)),rep('7',length(inds7)),
            rep('8',length(inds8)),rep('9',length(inds9)))
            
DElist = data.frame( DEgenes = DElist, labels = as.factor(labels))
write.csv(DElist, file = "~/Desktop/Eva/OnClusterPlots/DElist.csv",row.names=FALSE)

#Next step ...

xenopus.combined$celltype.stim <- paste(Idents(xenopus.combined), xenopus.combined$stim, sep = "_")
xenopus.combined$celltype <- Idents(xenopus.combined)
Idents(xenopus.combined) <- "celltype.stim"
NT.response0 <- FindMarkers(xenopus.combined, ident.1 = "0_SIGAH12", ident.2 = "0_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response0), file = "~/Desktop/Eva/OnClusterPlots/DElist_0.csv")
NT.response1 <- FindMarkers(xenopus.combined, ident.1 = "1_SIGAH12", ident.2 = "1_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response1), file = "~/Desktop/Eva/OnClusterPlots/DElist_1.csv")
NT.response2 <- FindMarkers(xenopus.combined, ident.1 = "2_SIGAH12", ident.2 = "2_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response2), file = "~/Desktop/Eva/OnClusterPlots/DElist_2.csv")
NT.response3 <- FindMarkers(xenopus.combined, ident.1 = "3_SIGAH12", ident.2 = "3_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response3), file = "~/Desktop/Eva/OnClusterPlots/DElist_3.csv")
NT.response4 <- FindMarkers(xenopus.combined, ident.1 = "4_SIGAH12", ident.2 = "4_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response4), file = "~/Desktop/Eva/OnClusterPlots/DElist_4.csv")
NT.response5 <- FindMarkers(xenopus.combined, ident.1 = "5_SIGAH12", ident.2 = "5_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response5), file = "~/Desktop/Eva/OnClusterPlots/DElist_5.csv")
NT.response6 <- FindMarkers(xenopus.combined, ident.1 = "6_SIGAH12", ident.2 = "6_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response6), file = "~/Desktop/Eva/OnClusterPlots/DElist_6.csv")
NT.response7 <- FindMarkers(xenopus.combined, ident.1 = "7_SIGAH12", ident.2 = "7_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response7), file = "~/Desktop/Eva/OnClusterPlots/DElist_7.csv")
NT.response8 <- FindMarkers(xenopus.combined, ident.1 = "8_SIGAH12", ident.2 = "8_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response8), file = "~/Desktop/Eva/OnClusterPlots/DElist_8.csv")
NT.response9 <- FindMarkers(xenopus.combined, ident.1 = "9_SIGAH12", ident.2 = "9_SIGAH5", verbose = FALSE)
write.csv(as.data.frame(NT.response9), file = "~/Desktop/Eva/OnClusterPlots/DElist_9.csv")
#NT.response10 <- FindMarkers(xenopus.combined, ident.1 = "10_SIGAH12", ident.2 = "10_SIGAH5", verbose = FALSE)
#write.csv(as.data.frame(NT.response10), file = "~/Desktop/Eva/DOnClusterPlots/Elist_10.csv")


#head(NT.response, n = 15)
for (i in 1:20){
FeaturePlot(xenopus.combined, features = rownames(NT.response0)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl0_DEplots_",i,".pdf",sep=""))
}

for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response1)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl1_DEplots_",i,".pdf",sep=""))
}

for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response2)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl2_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response3)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl3_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response4)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl4_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response5)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl5_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response6)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl6_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response7)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl7_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response8)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl8_DEplots_",i,".pdf",sep=""))
}
for (i in 1:20){
  FeaturePlot(xenopus.combined, features = rownames(NT.response9)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","Cl9_DEplots_",i,".pdf",sep=""))
}
#for (i in 1:20){
#  FeaturePlot(xenopus.combined, features = rownames(NT.response10)[i], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
#  ggsave(filename=paste("~/Desktop/Eva/ClusterPlots/OnClusterPlots/","Cl10_DEplots_",i,".pdf",sep=""))
#}



#FeaturePlot(xenopus.combined, features = OnMem[1,], split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
OnMem <- read.delim("/Users/christopherpenfold/Desktop/Eva/Motifs/Genelists_from_Eva_April2018/EndoNT_Filter_DonorsvsIVFdownFDR_DonsRPKM1_IVFvsNT_FDRFC3up__OnMem3FC.txt")

for (i in 1:dim(OnMem)[1]){
FeaturePlot(xenopus.combined, features = as.character(OnMem[i,]), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","OnMem_",as.character(OnMem[i,]),".pdf",sep=""))
}

OffMem <- read.delim("/Users/christopherpenfold/Desktop/Eva/Motifs/Genelists_from_Eva_April2018/EndoNT_Filter_DonorsvsIVFupFDR_IVFvsNT_FDRFC3down__OffMem3FC.txt")

for (i in 1:dim(OffMem)[1]){
  FeaturePlot(xenopus.combined, features = as.character(OffMem[i,]), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","OffMem_",as.character(OnMem[i,]),".pdf",sep=""))
}


for (i in 1:dim(OffMem)[1]){
  plots <- VlnPlot(xenopus.combined, features = as.character(OffMem[i,]), split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
  CombinePlots(plots = plots, ncol = 1)
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","OffMem_",as.character(OnMem[i,]),"_violin.pdf",sep=""))
}

for (i in 1:dim(OnMem)[1]){
  plots <- VlnPlot(xenopus.combined, features = as.character(OnMem[i,]), split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
  CombinePlots(plots = plots, ncol = 1)
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/","OnMem_",as.character(OnMem[i,]),"_violin.pdf",sep=""))
}

#Markers <- read.delim("/Users/christopherpenfold/Desktop/Eva/Markers.txt")

library(stringr)

for (i in 1:length(levels(Markers$State))){
  
  tryCatch({
  
  currstateind <- which(Markers$State==levels(Markers$State)[i])
  markerlist <- c(paste(Markers$Marker_genes[currstateind],".L",sep = ""),paste(Markers$Marker_genes[currstateind],".S",sep = ""))
  filename <- levels(Markers$State)[i]
  filename <- str_replace_all(filename, fixed(" "), "_")
  filename <- str_replace_all(filename, fixed("/"), "_")
  filename <- str_replace_all(filename, fixed("("), "_")
  filename <- str_replace_all(filename, fixed(")"), "_")
  sf <- max(c(1,length(markerlist)/6))
  FeaturePlot(xenopus.combined, features = markerlist, min.cutoff = "q9")
  ggsave(filename=paste("~/Desktop/Eva/OnClusterPlots/",filename,".pdf",sep=""), width = 13.7*sf, height = 4.86*sf,limitsize = FALSE)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

