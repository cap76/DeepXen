#Innstall reticulate (uncommet the line below to do so)
#install.packages("reticulate")
library(reticulate)
library(pheatmap)
library(ggplot2)

setwd("/Volumes/Overflow")
np <- import("numpy")

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("SomEndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("SomEndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("SomEndoFull/On_attributions_OnTP_perm.npy")

#Index for TP ON memory genes (to get get IDS etc)
onTP <- np$load("SomEndoFull/OnTP.npy")

#Get the gene IDS etc
Data <-read.table("SomEndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- as.data.frame(Data$V4[onTP])
colnames(genes) <- c("ID")

TFList <- read.table("/Users/christopherpenfold/Desktop/Eva/TFList.csv",sep=",", header = T, row.names=NULL)

AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)


#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
redblue1 <- colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

#Do k-means on average histone values and then look  at the overall profiles.
#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

colnames(mean1) <- labs
NCl <- 4
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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl_delta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl1_delta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl2_delta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl3_delta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl4_delta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl5_delta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl6_delta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl7_delta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl8_delta.pdf")


ONList1 <- unique(as.character(AllG[ which(clusters$cluster==1) ,2]))
ONList2 <- unique(as.character(AllG[ which(clusters$cluster==2) ,2]))
ONList3 <- unique(as.character(AllG[ which(clusters$cluster==3) ,2]))
ONList4 <- unique(as.character(AllG[ which(clusters$cluster==4) ,2]))
ONList5 <- unique(as.character(AllG[ which(clusters$cluster==5) ,2]))
ONList6 <- unique(as.character(AllG[ which(clusters$cluster==6) ,2]))
ONList7 <- unique(as.character(AllG[ which(clusters$cluster==7) ,2]))
ONList8 <- unique(as.character(AllG[ which(clusters$cluster==8) ,2]))

saveext <- "./"
write.csv(ONList1, file=paste(saveext,"/SomEndo_Heatmaps/ONList1.csv",sep=""))
write.csv(ONList2, file=paste(saveext,"/SomEndo_Heatmaps/ONList2.csv",sep=""))
write.csv(ONList3, file=paste(saveext,"/SomEndo_Heatmaps/ONList3.csv",sep=""))
write.csv(ONList4, file=paste(saveext,"/SomEndo_Heatmaps/ONList4.csv",sep=""))
write.csv(ONList5, file=paste(saveext,"/SomEndo_Heatmaps/ONList5.csv",sep=""))
write.csv(ONList6, file=paste(saveext,"/SomEndo_Heatmaps/ONList6.csv",sep=""))
write.csv(ONList7, file=paste(saveext,"/SomEndo_Heatmaps/ONList7.csv",sep=""))
write.csv(ONList8, file=paste(saveext,"/SomEndo_Heatmaps/ONList8.csv",sep=""))

lowerCI <- function(vec) {
  inds <- as.integer(length(vec)*0.05)
  vec <- sort(vec)[inds]
  return(vec)
}

upperCI <- function(vec) {
  inds <- as.integer(length(vec)*0.95)
  vec <- sort(vec)[inds]
  return(vec)
}

#Now plot an example gene ...
H3K4me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3 <- H3K4me3[onTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me3_1 <- H3K4me3[ind,]
  M <- colMeans(H3K4me3_1)
  V <- apply(H3K4me3_1, c(2), var) 
  L <- apply(H3K4me3_1, c(2), lowerCI) 
  U <- apply(H3K4me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K4me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H2Aub <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H2A.ub.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H2Aub <- H2Aub[onTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H2Aub_1 <- H2Aub[ind,]
  M <- colMeans(H2Aub_1)
  V <- apply(H2Aub_1, c(2), var) 
  L <- apply(H2Aub_1, c(2), lowerCI) 
  U <- apply(H2Aub_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H2Aub_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K4me1 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K4me1.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me1 <- H3K4me1[onTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me1_1 <- H3K4me1[ind,]
  M <- colMeans(H3K4me1_1)
  V <- apply(H3K4me1_1, c(2), var) 
  L <- apply(H3K4me1_1, c(2), lowerCI) 
  U <- apply(H3K4me1_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K4me1_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K9me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K9me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K9me3 <- H3K9me3[onTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K9me3_1 <- H3K9me3[ind,]
  M <- colMeans(H3K9me3_1)
  V <- apply(H3K9me3_1, c(2), var) 
  L <- apply(H3K9me3_1, c(2), lowerCI) 
  U <- apply(H3K9me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K9me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K27ac <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K27ac.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27ac <- H3K27ac[onTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27ac_1 <- H3K27ac[ind,]
  M <- colMeans(H3K27ac_1)
  V <- apply(H3K27ac_1, c(2), var) 
  L <- apply(H3K27ac_1, c(2), lowerCI) 
  U <- apply(H3K27ac_1, c(2), upperCI)   
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K27ac_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K27me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K27me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27me3 <- H3K27me3[onTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27me3_1 <- H3K27me3[ind,]
  M <- colMeans(H3K27me3_1)
  V <- apply(H3K27me3_1, c(2), var) 
  L <- apply(H3K27me3_1, c(2), lowerCI) 
  U <- apply(H3K27me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K27me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K36me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K36me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K36me3 <- H3K36me3[onTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K36me3_1 <- H3K36me3[ind,]
  M <- colMeans(H3K36me3_1)
  V <- apply(H3K36me3_1, c(2), var) 
  L <- apply(H3K36me3_1, c(2), lowerCI) 
  U <- apply(H3K36me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K36me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K79me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K79me3.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K79me3 <- H3K79me3[onTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K79me3_1 <- H3K79me3[ind,]
  M <- colMeans(H3K79me3_1)
  V <- apply(H3K79me3_1, c(2), var) 
  L <- apply(H3K79me3_1, c(2), lowerCI) 
  U <- apply(H3K79me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K79me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}


for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl_nodelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl1_nodelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl2_nodelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl3_nodelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl4_nodelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl5_nodelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl6_nodelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl7_nodelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/ONRDPCl8_nodelta.pdf")





#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("SomEndoFull/Up_attributions_OffTP.npy")
mat1 <- np$load("SomEndoFull/Off_attributions_OffTP.npy")
mat2 <- np$load("SomEndoFull/On_attributions_OffTP_perm.npy")

#Index for TP ON memory genes (to get get IDS etc)
offTP <- np$load("SomEndoFull/OffTP.npy")

#Get the gene IDS etc
Data <-read.table("SomEndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- as.data.frame(Data$V4[offTP])
colnames(genes) <- c("ID")

TFList <- read.table("/Users/christopherpenfold/Desktop/Eva/TFList.csv",sep=",", header = T, row.names=NULL)

AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)

labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
redblue1 <- colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))



#p(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFF.pdf")

mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
#mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

colnames(mean1) <- labs
NCl <- 4
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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl_delta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl1_delta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFNRUPCl2_delta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl3_delta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl4_delta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl5_delta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl6_delta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl7_delta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRDUCl8_delta.pdf")


OFFList1 <- unique(as.character(AllG[ which(clusters$cluster==1) ,2]))
OFFList2 <- unique(as.character(AllG[ which(clusters$cluster==2) ,2]))
OFFList3 <- unique(as.character(AllG[ which(clusters$cluster==3) ,2]))
OFFList4 <- unique(as.character(AllG[ which(clusters$cluster==4) ,2]))
OFFList5 <- unique(as.character(AllG[ which(clusters$cluster==5) ,2]))
OFFList6 <- unique(as.character(AllG[ which(clusters$cluster==6) ,2]))
OFFList7 <- unique(as.character(AllG[ which(clusters$cluster==7) ,2]))
OFFList8 <- unique(as.character(AllG[ which(clusters$cluster==8) ,2]))

saveext <- "./"
write.csv(OFFList1, file=paste(saveext,"/SomEndo_Heatmaps/OFFList1.csv",sep=""))
write.csv(OFFList2, file=paste(saveext,"/SomEndo_Heatmaps/OFFList2.csv",sep=""))
write.csv(OFFList3, file=paste(saveext,"/SomEndo_Heatmaps/OFFList3.csv",sep=""))
write.csv(OFFList4, file=paste(saveext,"/SomEndo_Heatmaps/OFFList4.csv",sep=""))
write.csv(OFFList5, file=paste(saveext,"/SomEndo_Heatmaps/OFFList5.csv",sep=""))
write.csv(OFFList6, file=paste(saveext,"/SomEndo_Heatmaps/OFFList6.csv",sep=""))
write.csv(OFFList7, file=paste(saveext,"/SomEndo_Heatmaps/OFFList7.csv",sep=""))
write.csv(OFFList8, file=paste(saveext,"/SomEndo_Heatmaps/OFFList8.csv",sep=""))




#Now plot an example gene ...
H3K4me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3 <- H3K4me3[offTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me3_1 <- H3K4me3[ind,]
  M <- colMeans(H3K4me3_1)
  V <- apply(H3K4me3_1, c(2), var) 
  L <- apply(H3K4me3_1, c(2), lowerCI) 
  U <- apply(H3K4me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K4me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H2Aub <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H2A.ub.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H2Aub <- H2Aub[offTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H2Aub_1 <- H2Aub[ind,]
  M <- colMeans(H2Aub_1)
  V <- apply(H2Aub_1, c(2), var) 
  L <- apply(H2Aub_1, c(2), lowerCI) 
  U <- apply(H2Aub_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H2Aub_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K4me1 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K4me1.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me1 <- H3K4me1[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me1_1 <- H3K4me1[ind,]
  M <- colMeans(H3K4me1_1)
  V <- apply(H3K4me1_1, c(2), var) 
  L <- apply(H3K4me1_1, c(2), lowerCI) 
  U <- apply(H3K4me1_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K4me1_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K9me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K9me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K9me3 <- H3K9me3[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K9me3_1 <- H3K9me3[ind,]
  M <- colMeans(H3K9me3_1)
  V <- apply(H3K9me3_1, c(2), var) 
  L <- apply(H3K9me3_1, c(2), lowerCI) 
  U <- apply(H3K9me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K9me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K27ac <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K27ac.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27ac <- H3K27ac[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27ac_1 <- H3K27ac[ind,]
  M <- colMeans(H3K27ac_1)
  V <- apply(H3K27ac_1, c(2), var) 
  L <- apply(H3K27ac_1, c(2), lowerCI) 
  U <- apply(H3K27ac_1, c(2), upperCI)   
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K27ac_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K27me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K27me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27me3 <- H3K27me3[offTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27me3_1 <- H3K27me3[ind,]
  M <- colMeans(H3K27me3_1)
  V <- apply(H3K27me3_1, c(2), var) 
  L <- apply(H3K27me3_1, c(2), lowerCI) 
  U <- apply(H3K27me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K27me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K36me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K36me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K36me3 <- H3K36me3[offTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K36me3_1 <- H3K36me3[ind,]
  M <- colMeans(H3K36me3_1)
  V <- apply(H3K36me3_1, c(2), var) 
  L <- apply(H3K36me3_1, c(2), lowerCI) 
  U <- apply(H3K36me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K36me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

H3K79me3 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K79me3.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K79me3 <- H3K79me3[offTP,]

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K79me3_1 <- H3K79me3[ind,]
  M <- colMeans(H3K79me3_1)
  V <- apply(H3K79me3_1, c(2), var) 
  L <- apply(H3K79me3_1, c(2), lowerCI) 
  U <- apply(H3K79me3_1, c(2), upperCI)   
  
  D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("SomEndo_Heatmaps/H3K79me3_OFFPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}


SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(NCl,dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl_nodelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl1_nodelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFNRUPCl2_nodelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl3_nodelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl4_nodelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl5_nodelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl6_nodelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRUPCl7_nodelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomEndo_Heatmaps/OFFRDUCl8_nodelta.pdf")



sadsadsadadadsaad




#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.#
#mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "SomeEndo_OFF.pdf")

mean1 = apply(mat1-mat, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

#colnames(mean1) <- labs
NCl <- 4
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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl_Delta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl1_Delta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFNRUPCl2_Delta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl3_Delta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl4_Delta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl5_Delta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl6_Delta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl7_Delta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRDUCl8_Delta.pdf")

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl1.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFNRUPCl2.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl3.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl4.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl5.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl6.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRUPCl7.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "SomeEndo_OFFRDUCl8.pdf")




OFFList1 <- unique(as.character(AllG[ which(clusters$cluster==1) ,2]))
OFFList2 <- unique(as.character(AllG[ which(clusters$cluster==2) ,2]))
OFFList3 <- unique(as.character(AllG[ which(clusters$cluster==3) ,2]))
OFFList4 <- unique(as.character(AllG[ which(clusters$cluster==4) ,2]))
OFFList5 <- unique(as.character(AllG[ which(clusters$cluster==5) ,2]))
OFFList6 <- unique(as.character(AllG[ which(clusters$cluster==6) ,2]))
OFFList7 <- unique(as.character(AllG[ which(clusters$cluster==7) ,2]))
OFFList8 <- unique(as.character(AllG[ which(clusters$cluster==8) ,2]))

saveext <- "./"
write.csv(OFFList1, file=paste(saveext,"/SomEndo_OFFList1.csv",sep=""))
write.csv(OFFList2, file=paste(saveext,"/SomEndo_OFFList2.csv",sep=""))
write.csv(OFFList3, file=paste(saveext,"/SomEndo_OFFList3.csv",sep=""))
write.csv(OFFList4, file=paste(saveext,"/SomEndo_OFFList4.csv",sep=""))
write.csv(OFFList5, file=paste(saveext,"/SomEndo_OFFList5.csv",sep=""))
write.csv(OFFList6, file=paste(saveext,"/SomEndo_OFFList6.csv",sep=""))
write.csv(OFFList7, file=paste(saveext,"/SomEndo_OFFList7.csv",sep=""))
write.csv(OFFList8, file=paste(saveext,"/SomEndo_OFFList8.csv",sep=""))


#Now plot an example gene ...
H3K4me3_1 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Ecto_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3_1 <- H3K4me3_1[offTP,]

H3K4me3_2 <- read.table("/Volumes/Overflow/HistoneMatrix/SomEndo/Endo_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3_2 <- H3K4me3_2[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me3_1b <- H3K4me3_1[ind,]
  M <- colMeans(H3K4me3_1b)
  V <- apply(H3K4me3_1b, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)

  H3K4me3_2b <- H3K4me3_2[ind,]
  M <- colMeans(H3K4me3_2b)
  V <- apply(H3K4me3_2b, c(2), var) 
  D2 <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  
  
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K4me3_SomeEndo_OFFCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K36me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K36me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K36me3 <- H3K36me3[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K36me3_1 <- H3K36me3[ind,]
  M <- colMeans(H3K36me3_1)
  V <- apply(H3K4me3_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K36me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K27ac <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K27ac.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27ac <- H3K27ac[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27ac_1 <- H3K27ac[ind,]
  M <- colMeans(H3K27ac_1)
  V <- apply(H3K27ac_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K27ac_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K79me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K79me3.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K79me3 <- H3K79me3[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K79me3_1 <- H3K79me3[ind,]
  M <- colMeans(H3K79me3_1)
  V <- apply(H3K79me3_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K79me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H2Aub <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H2A.ub.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H2Aub <- H2Aub[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H2Aub_1 <- H2Aub[ind,]
  M <- colMeans(H2Aub_1)
  V <- apply(H2Aub_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H2Aub_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K9me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K9me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K9me3 <- H3K9me3[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K9me3_1 <- H3K9me3[ind,]
  M <- colMeans(H3K9me3_1)
  V <- apply(H3K9me3_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K9me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}


#Now plot an example gene ...
H3K4me1 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K4me1.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me1 <- H3K4me1[offTP,]
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me1_1 <- H3K4me1[ind,]
  M <- colMeans(H3K4me1_1)
  V <- apply(H3K4me1_1, c(2), var) 
  D <- data.frame(X=seq(1,length(M),1),Y=M,Z=V)
  p <- ggplot(D, aes(X, Y))+geom_line(aes(x=X, y=Y), size = 1)+geom_ribbon(aes(x=X, ymin=M-V, ymax=M+V), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("H3K4me1_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}






sdsadasdsda

#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ON.pdf")

mean2 = apply(mat1[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat1[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat1[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ONTF.pdf")

mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RD.pdf")

mean2 = apply(mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RDTF.pdf")


mean2 = apply(mat1-mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONmiusRD.pdf")

mean2 = apply(mat1-mat2, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat2, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONminusPerm.pdf")



mean2 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONmiusRDTF.pdf")

mean2 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONminusPermTF.pdf")









#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("SomEndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("SomEndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("SomEndoFull/On_attributions_OnTP_perm.npy")

#Index for TP ON memory genes (to get get IDS etc)
offTP <- np$load("SomEndoFull/OffTP.npy")

#Get the gene IDS etc
Data <-read.table("SomEndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- as.data.frame(Data$V4[offTP])
colnames(genes) <- c("ID")

AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)

#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
redblue1 <- colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))


#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFF.pdf")

mean2 = apply(mat1[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat1[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat1[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFFTF.pdf")

mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RU.pdf")

mean2 = apply(mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RUTF.pdf")


mean2 = apply(mat1-mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFmiusRU.pdf")

mean2 = apply(mat1-mat2, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat2, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFminusPerm.pdf")



mean2 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFmiusRUTF.pdf")

mean2 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
mean3 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFminusPermTF.pdf")

