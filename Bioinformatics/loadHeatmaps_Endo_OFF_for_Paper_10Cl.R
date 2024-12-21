#Innstall reticulate (uncommet the line below to do so)
#install.packages("reticulate")
#library(reticulate)
#library(pheatmap)

library(reticulate)
use_python("/usr/bin/python3", required = T)
sys = import("sys")

setwd("/Users/christopherpenfold/Desktop/Code/DeepXen/")
#setwd("/Volumes/Overflow")
np <- import("numpy")

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoFull/Up_attributions_OffTP.npy")
mat1 <- np$load("EndoFull/Off_attributions_OffTP.npy")
mat2 <- np$load("EndoFull/Off_attributions_OffTP_perm.npy")

#mat1 <- np$load("EndoFull/Down_attributions_OnTP_1.npy")


#Index for TP ON memory genes (to get get IDS etc)
offTP <- np$load("EndoFull/OffTP.npy")



#Get the gene IDS etc
Data <-read.table("EndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- Data$V4[offTP]
#colnames(genes) <- c("ID")
#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

genes


#Do k-means on average histone values and then look  at the overall profiles.
#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

#"#f0f9ff","#f7f0f0"
redblue1 <- colorRampPalette(c("#00B0F0",
                               "#DFF5FD",
                               "#FFFFFF",
                               "#FFEBEA",
                               "#FF0B07"))

#redblue1 <- colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

#"#f7f0ff"
#"#FFE1E0"
#"#FFF5F5"
#"#FFEBEA"
#"#FFCDCC"
#"#FFD7D6"
#"#FFC2C1"
#"#FFB8B7"

#"#f0f9ff"
#"#F4FCFE"
#"#EAF8FE"
#"#DFF5FD"
#"#D5F2FC"
#"#CAEFFC"
##C0EBFB

colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

pheatmap(t(SummaryMat1), color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_delta_defaultscale.pdf")



#"#f7f0ff"
#"#FFE1E0"
#"#FFF5F5"
#"#FFEBEA"
#"#FFCDCC"
#"#FFD7D6"
#"#FFC2C1"
#"#FFB8B7"

#"#f0f9ff"
#"#F4FCFE"
#"#EAF8FE"
#"#DFF5FD"
#"#D5F2FC"
#"#CAEFFC"
##C0EBFB
redblue1 <- colorRampPalette(c("#00B0F0",
                               "#D5F2FC",
                               "#FFFFFF",
                               "#FFB8B7",
                               "#FF0B07"))

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_delta_scalein0p2.pdf")

mat_breaks <- seq(-0.1, 0.1, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_delta_scalein0p1.pdf")

mat_breaks <- seq(-0.2, 0.6, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_delta_scalein-0p20p6.pdf")

#mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)
#setwd("/Users/christopherpenfold/Desktop/Code/DeepXen/")


mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1, c(2,3), sum) / dim(mat1)[1]

mean3 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean4 = apply(mat, c(2,3), sum) / dim(mat)[1]

mean5 = apply(mat2, c(1,3), sum) / dim(mat2)[2]
mean6 = apply(mat2, c(2,3), sum) / dim(mat2)[1]

SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat2 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat3 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat2 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat3 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
pVal1 <- array(0,c(NCl,dim(mat)[3]))
pVal2 <- array(0,c(NCl,dim(mat)[3]))

#mean1[ind,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  uniqueMat2[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(ind)
  uniqueMat3[i,,] <- apply(mat2[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- colSums(mean1[ind,])/length(ind)
  SummaryMat2[i,] <- colSums(mean3[ind,])/length(ind)
  SummaryMat3[i,] <- colSums(mean5[ind,])/length(ind)
  
  for (j in 1:22) {
    p1 <- t.test(mean1[ind,j],mean3[ind,])
    p2 <- t.test(mean1[ind,j],mean5[ind,])
    pVal1[i,j] <- p1$p.value
    pVal2[i,j] <- p2$p.value
  }
}
SummaryMat <- rbind(SummaryMat1,SummaryMat2,SummaryMat3)
colnames(SummaryMat) <- labs
rownames(SummaryMat) <- c(paste("Cl",rep(1:NCl, each=1)), paste("ClR",rep(1:NCl, each=1)), paste("ClP",rep(1:NCl, each=1)) )
pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_default.pdf")
colnames(pVal1) <- labs
colnames(pVal2) <- labs

rownames(pVal1) <- c(paste("Cl",rep(1:NCl, each=1)) )
rownames(pVal2) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(0, 0.05, length.out = 20)
pheatmap(t(pVal1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_default_pvaluevsmemory.pdf")
pheatmap(t(pVal2),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_default_pvaluevsshuffle.pdf")


mat_breaks <- seq(-0.0004, 0.0009, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20), breaks=mat_breaks,gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_0p0004_0p0009.pdf")

mat_breaks <- seq(-0.0002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20), breaks=mat_breaks,gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_0p0002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0004, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20), breaks=mat_breaks,gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_0p0002_0p0004.pdf")


mat_breaks <- seq(-0.0002, 0.0003, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20), breaks=mat_breaks,gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_0p0002_0p0003.pdf")

mat_breaks <- seq(-0.0002, 0.0002, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20), breaks=mat_breaks,gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_0p0002_0p0002.pdf")

#write.table(data.frame(genes,clusters$cluster),file="~/Desktop/Clusters_Endo_Off.csv")

#
mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]
colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/OFFRUPCl_delta_10Cl_scalein0p2.pdf")



#No On vs Downs
redblue1 <- colorRampPalette(c("#00B0F0",
                               "#D5F2FC",
                               "#FFFFFF",
                               "#FFB8B7",
                               "#FF0B07"))
#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("EndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("EndoFull/On_attributions_OnTP_perm.npy")

#Index for TP ON memory genes (to get get IDS etc)
onTP <- np$load("EndoFull/OnTP.npy")

#Get the gene IDS etc
Data <-read.table("EndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
#genes <- as.data.frame(Data$V4[onTP])
genes <- Data$V4[onTP]

#colnames(genes) <- c("ID")

#AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)

#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")


#Example to print out the heatmap and save it
#mat_breaks <- seq(0, 2, length.out = 20) #Set the scale bar manually rather than automatically
#pheatmap(mat,color =  redblue1(20),  border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = "HMexample.pdf")

#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.

mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]



#mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ON.pdf")


colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

pheatmap(t(SummaryMat1), color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_delta_defaultscale.pdf")


mat_breaks <- seq(-0.1, 0.4, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/ONRDCl_delta_scalein0p1_0p4.pdf")


mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/ONRDCl_delta_scalein0p2.pdf")


mat_breaks <- seq(-0.1, 0.1, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/ONRDCl_delta_scalein0p1.pdf")


mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1, c(2,3), sum) / dim(mat1)[1]

mean3 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean4 = apply(mat, c(2,3), sum) / dim(mat)[1]

mean5 = apply(mat2, c(1,3), sum) / dim(mat2)[2]
mean6 = apply(mat2, c(2,3), sum) / dim(mat2)[1]

SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat2 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat3 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat2 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat3 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
pVal1 <- array(0,c(NCl,dim(mat)[3]))
pVal2 <- array(0,c(NCl,dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  uniqueMat2[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(ind)
  uniqueMat3[i,,] <- apply(mat2[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- colSums(mean1[ind,])/length(ind)
  SummaryMat2[i,] <- colSums(mean3[ind,])/length(ind)
  SummaryMat3[i,] <- colSums(mean5[ind,])/length(ind)
  for (j in 1:22) {
    p1 <- t.test(mean1[ind,j],mean3[ind,])
    p2 <- t.test(mean1[ind,j],mean5[ind,])
    pVal1[i,j] <- p1$p.value
    pVal2[i,j] <- p2$p.value
  }
}
SummaryMat <- rbind(SummaryMat1,SummaryMat2,SummaryMat3)
colnames(SummaryMat) <- labs
rownames(SummaryMat) <- c(paste("Cl",rep(1:NCl, each=1)), paste("ClR",rep(1:NCl, each=1)), paste("ClP",rep(1:NCl, each=1)) )
pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_default.pdf")

mat_breaks <- seq(-0.0002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_scalein0p0002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0004, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_scalein0p0002_0p0004.pdf")



mat_breaks <- seq(-0.0001, 0.0002, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_scalein0p0001_0p0002.pdf")

write.table(data.frame(genes,clusters$cluster),file="~/Desktop/Clusters_Endo_On.csv")

#
mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]
colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_delta_10Cl_scalein0p2.pdf")

rownames(pVal1) <- c(paste("Cl",rep(1:NCl, each=1)) )
rownames(pVal2) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(0, 0.05, length.out = 20)
pheatmap(t(pVal1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_default_pvaluevsmemory.pdf")
pheatmap(t(pVal2),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/ONRDPCl_default_pvaluevsshuffle.pdf")

#Now somite endo ON

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("SomEndoFull/Up_attributions_OffTP.npy")
mat1 <- np$load("SomEndoFull/Off_attributions_OffTP.npy")
mat2 <- np$load("SomEndoFull/Off_attributions_OffTP_perm.npy")

offTP <- np$load("SomEndoFull/OffTP.npy")

#Get the gene IDS etc
Data <-read.table("SomEndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- Data$V4[offTP]
#colnames(genes) <- c("ID")
#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")
#labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

pheatmap(t(SummaryMat1), color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_OFFRUPCl_delta_defaultscale.pdf")


mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_OFFRUCl_delta_scalein0p2.pdf")


mat_breaks <- seq(-0.1, 0.1, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_OFFRUCl_delta_scalein0p1.pdf")



mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1, c(2,3), sum) / dim(mat1)[1]

mean3 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean4 = apply(mat, c(2,3), sum) / dim(mat)[1]

mean5 = apply(mat2, c(1,3), sum) / dim(mat2)[2]
mean6 = apply(mat2, c(2,3), sum) / dim(mat2)[1]

SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat2 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat3 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat2 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat3 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
pVal1 <- array(0,c(NCl,dim(mat)[3]))
pVal2 <- array(0,c(NCl,dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  uniqueMat2[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(ind)
  uniqueMat3[i,,] <- apply(mat2[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- colSums(mean1[ind,])/length(ind)
  SummaryMat2[i,] <- colSums(mean3[ind,])/length(ind)
  SummaryMat3[i,] <- colSums(mean5[ind,])/length(ind)
  for (j in 1:22) {
    p1 <- t.test(mean1[ind,j],mean3[ind,])
    p2 <- t.test(mean1[ind,j],mean5[ind,])
    pVal1[i,j] <- p1$p.value
    pVal2[i,j] <- p2$p.value
  }
}
SummaryMat <- rbind(SummaryMat1,SummaryMat2,SummaryMat3)
colnames(SummaryMat) <- labs
rownames(SummaryMat) <- c(paste("Cl",rep(1:NCl, each=1)), paste("ClR",rep(1:NCl, each=1)), paste("ClP",rep(1:NCl, each=1)) )
pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_OFFRUPCl_default.pdf")

rownames(pVal1) <- c(paste("Cl",rep(1:NCl, each=1)) )
rownames(pVal2) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(0, 0.05, length.out = 20)
pheatmap(t(pVal1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomeEndoOFFRUPCl_default_pvaluevsmemory.pdf")
pheatmap(t(pVal2),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndoOFFRUPCl_default_pvaluevsshuffle.pdf")


mat_breaks <- seq(-0.002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/omEndo_OFFRUPCl_scalein0p002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/omEndo_OFFRUPCl_scalein0p0002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0004, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/omEndo_OFFRUPCl_scalein0p0002_0p0004.pdf")



mat_breaks <- seq(-0.0002, 0.0002, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/omEndo_OFFRUPCl_scalein0p0002_0p0002.pdf")




mat <- np$load("SomEndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("SomEndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("SomEndoFull/On_attributions_OnTP_perm.npy")



onTP <- np$load("SomEndoFull/OnTP.npy")

#Get the gene IDS etc
Data <-read.table("SomEndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- Data$V4[onTP]
#colnames(genes) <- c("ID")
#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")
#labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

colnames(mean1) <- labs
NCl <- 10
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

pheatmap(t(SummaryMat1), color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDPCl_delta_defaultscale.pdf")


mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDCl_delta_scalein0p2.pdf")


mat_breaks <- seq(-0.1, 0.1, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDCl_delta_scalein0p1.pdf")


mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1, c(2,3), sum) / dim(mat1)[1]

mean3 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean4 = apply(mat, c(2,3), sum) / dim(mat)[1]

mean5 = apply(mat2, c(1,3), sum) / dim(mat2)[2]
mean6 = apply(mat2, c(2,3), sum) / dim(mat2)[1]

SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat2 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat3 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat2 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat3 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
pVal1 <- array(0,c(NCl,dim(mat)[3]))
pVal2 <- array(0,c(NCl,dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  uniqueMat2[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(ind)
  uniqueMat3[i,,] <- apply(mat2[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- colSums(mean1[ind,])/length(ind)
  SummaryMat2[i,] <- colSums(mean3[ind,])/length(ind)
  SummaryMat3[i,] <- colSums(mean5[ind,])/length(ind)
  for (j in 1:22) {
    p1 <- t.test(mean1[ind,j],mean3[ind,])
    p2 <- t.test(mean1[ind,j],mean5[ind,])
    pVal1[i,j] <- p1$p.value
    pVal2[i,j] <- p2$p.value
  }
}
SummaryMat <- rbind(SummaryMat1,SummaryMat2,SummaryMat3)
colnames(SummaryMat) <- labs
rownames(SummaryMat) <- c(paste("Cl",rep(1:NCl, each=1)), paste("ClR",rep(1:NCl, each=1)), paste("ClP",rep(1:NCl, each=1)) )
pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDPCl_default.pdf")


rownames(pVal1) <- c(paste("Cl",rep(1:NCl, each=1)) )
rownames(pVal2) <- c(paste("Cl",rep(1:NCl, each=1)) )
mat_breaks <- seq(0, 0.05, length.out = 20)
pheatmap(t(pVal1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndoONRDPCl_default_pvaluevsmemory.pdf")
pheatmap(t(pVal2),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomeEndoONRDPCl_default_pvaluevsshuffle.pdf")


mat_breaks <- seq(-0.002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDPCl_scalein0p002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0006, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDPCl_scalein0p0002_0p0006.pdf")


mat_breaks <- seq(-0.0002, 0.0004, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_ONRDPCl_scalein0p0002_0p0004.pdf")



mat_breaks <- seq(-0.0002, 0.0002, length.out = 20)
pheatmap(t(SummaryMat),color =  redblue1(20),breaks = mat_breaks,  gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/omEndo_ONRDPCl_scalein0p0002_0p0002.pdf")



#mat_breaks <- seq(-0.0004, 0.0009, length.out = 20)
#pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "~/Desktop/SomEndo_OFFRUPCl_default.pdf")







