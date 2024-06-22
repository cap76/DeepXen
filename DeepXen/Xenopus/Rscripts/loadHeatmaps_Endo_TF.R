#Innstall reticulate (uncommet the line below to do so)
#install.packages("reticulate")
library(reticulate)
library(pheatmap)

setwd("/Volumes/Overflow")
np <- import("numpy")

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("EndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("EndoFull/Off_attributions_OnTP.npy")
mat3 <- np$load("EndoFull/Up_attributions_OnTP.npy")
mat4 <- np$load("EndoFull/Other_attributions_OnTP.npy")
mat5 <- np$load("EndoFull/On_attributions_OnTP_perm.npy")

#Now plot an example gene ...
#H3K4me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Ecto_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)

#Index for TP ON memory genes (to get get IDS etc)
onTP <- np$load("EndoFull/OnTP.npy")

#Compare onTP with DownTP?
rdTP <- np$load("EndoFull/RDTP.npy")
ruTP <- np$load("EndoFull/RUTP.npy")



 #onTP <- np$load("EndoFull/RDTP.npy")


TFList <- read.table("/Users/christopherpenfold/Desktop/Eva/TFList.csv",sep=",", header = T, row.names=NULL)

#Get the gene IDS etc
Data <-read.table("EndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- as.data.frame(Data$V4[onTP])
colnames(genes) <- c("ID")

genes2 <- as.data.frame(Data$V4[rdTP])
colnames(genes2) <- c("ID")


AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)
AllG2 <- merge(x = genes2, y = TFList, by = "ID", all.x = TRUE)


#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
redblue1 <- colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

#Plot example
mat_breaks <- seq(-0.008, 0.008, length.out = 20)
#colnames(mean3) <- labs
inds <- which(AllG[onTP,][,2]=="foxc2.L")[1]

x <- t(mat[inds,,])
rownames(x) <- labs
x1 <- t(mat1[inds,,])
rownames(x1) <- labs
x2 <- t(mat2[inds,,])
rownames(x2) <- labs
x3 <- t(mat3[inds,,])
rownames(x3) <- labs
x4 <- t(mat4[inds,,])
rownames(x4) <- labs
x5 <- t(mat5[inds,,])
rownames(x5) <- labs
pheatmap(x,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_down_FOXC2.pdf")
pheatmap(x1,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_on_FOXC2.pdf")
pheatmap(x2,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks,  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_off_FOXC2.pdf")
pheatmap(x3,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_up_FOXC2.pdf")
pheatmap(x4,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_other_FOXC2.pdf")
pheatmap(x5,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_perm_FOXC2.pdf")

#colnames(mean3) <- labs
x <- t(mat[15,,])
rownames(x) <- labs
x1 <- t(mat1[15,,])
rownames(x1) <- labs
x2 <- t(mat2[15,,])
rownames(x2) <- labs
x3 <- t(mat3[15,,])
rownames(x3) <- labs
x4 <- t(mat4[15,,])
rownames(x4) <- labs
x5 <- t(mat5[15,,])
rownames(x5) <- labs
pheatmap(x,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_down_wnt5a.pdf")
pheatmap(x1,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_on_wnt5a.pdf")
pheatmap(x2,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks,  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_off_wnt5a.pdf")
pheatmap(x3,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_up_wnt5a.pdf")
pheatmap(x4,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_other_wnt5a.pdf")
pheatmap(x5,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_perm_wnt5a.pdf")

#colnames(mean3) <- labs
inds <- which(AllG[onTP,][,2]=="kmt2a.L")[1]
x <- t(mat[inds,,])
rownames(x) <- labs
x1 <- t(mat1[inds,,])
rownames(x1) <- labs
x2 <- t(mat2[inds,,])
rownames(x2) <- labs
x3 <- t(mat3[inds,,])
rownames(x3) <- labs
x4 <- t(mat4[inds,,])
rownames(x4) <- labs
x5 <- t(mat5[inds,,])
rownames(x5) <- labs
pheatmap(x,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_down_kmt2a.pdf")
pheatmap(x1,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_on_kmt2a.pdf")
pheatmap(x2,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks,  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_off_kmt2a.pdf")
pheatmap(x3,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_up_kmt2a.pdf")
pheatmap(x4,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_other_kmt2a.pdf")
pheatmap(x5,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("..."), filename = "ONTP_perm_kmt2a.pdf")

#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
#mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ON.pdf")


library("grid")
library("gridExtra")

#H3 <- t(colSums(H3K4me3[onTP,]))/length(which(onTP=="TRUE"))

#ggplot(dat, aes(x=pos,y=value, colour=type)) +
#  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=type), alpha=0.3) +
#  theme_bw()

#y <- pheatmap(H3,color =  redblue1(20),  cluster_rows=FALSE,cluster_cols=FALSE, show_colnames = F, show_rownames = F,main = c("Average over all genes"))
#y <- ggplot(data=data.frame(x=seq(1,600,by=1),y=t(H3) ), aes(x=x, y=y)) + geom_line(color="blue")+ theme_classic() + scale_x_continuous(limits = c(0,600), expand = c(0, 0)) + ylab(NULL)
  
#x <- pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE)
#plot_list=list()
#plot_list[[1]] = y
#plot_list[[2]] = x[[4]]

#g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=1), widths = c(2), heights=c(1, 10))


#mean2 = apply(mat1[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat1[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat1[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ONTF.pdf")

#mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
#mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RD.pdf")

#mean2 = apply(mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RDTF.pdf")


#mean2 = apply(mat1-mat, c(1,3), sum) / dim(mat)[2]
#mean3 = apply(mat1-mat, c(2,3), sum) / dim(mat)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONmiusRD.pdf")

#mean2 = apply(mat1-mat2, c(1,3), sum) / dim(mat)[2]
##colnames(mean3) <- labs
#mean3 = apply(mat1-mat2, c(2,3), sum) / dim(mat)[1]
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONminusPerm.pdf")


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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl_delta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl1_delta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl2_delta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl3_delta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl4_delta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl5_delta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl6_delta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl7_delta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl8_delta.pdf")

#apply(as.matrix(df), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

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
H3K4me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3comp <- H3K4me3[rdTP,]
H3K4me3 <- H3K4me3[onTP,]
M1 <- colMeans(H3K4me3comp)
for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me3_1 <- H3K4me3[ind,]
  #M <- colMeans(H3K4me3_1)
  #V <- apply(H3K4me3_1, c(2), var) 
  #L <- apply(H3K4me3_1, c(2), lowerCI) 
  #U <- apply(H3K4me3_1, c(2), upperCI)   
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  
  M <- colMeans(H3K4me3_1)
  #V1 <- apply(H3K4me3_1, c(2), var) 
  #M2 <- colMeans(H3K4me3comp)
  #V2 <- apply(H3K4me3_1, c(2), var)   
  #L <- apply(H3K4me3_1, c(2), lowerCI) 
  #U <- apply(H3K4me3_1, c(2), upperCI)   
  
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                       ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
#    geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  
  
  ggsave(filename=paste("Endo_Heatmaps/H3K4me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K36me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K36me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K36me3comp <- H3K36me3[rdTP,]
H3K36me3 <- H3K36me3[onTP,]
M1 <- colMeans(H3K36me3comp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K36me3_1 <- H3K36me3[ind,]
  M <- colMeans(H3K36me3_1)
  #V <- apply(H3K36me3_1, c(2), var) 
  #L <- apply(H3K36me3_1, c(2), lowerCI) 
  #U <- apply(H3K36me3_1, c(2), upperCI)   
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  ggsave(filename=paste("Endo_Heatmaps/H3K36me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K27ac <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K27ac.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27accomp <- H3K27ac[rdTP,]
H3K27ac <- H3K27ac[onTP,]
M1 <- colMeans(H3K27accomp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27ac_1 <- H3K27ac[ind,]
  M <- colMeans(H3K27ac_1)
  #V <- apply(H3K27ac_1, c(2), var) 
  #L <- apply(H3K27ac_1, c(2), lowerCI) 
  #U <- apply(H3K27ac_1, c(2), upperCI)   
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K27ac_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K79me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K79me3.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K79me3comp <- H3K79me3[rdTP,]
H3K79me3 <- H3K79me3[onTP,]
M1 <- colMeans(H3K79me3comp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K79me3_1 <- H3K79me3[ind,]
  M <- colMeans(H3K79me3_1)
  #V <- apply(H3K79me3_1, c(2), var) 
  #L <- apply(H3K79me3_1, c(2), lowerCI) 
  #U <- apply(H3K79me3_1, c(2), upperCI)   
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K79me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H2Aub <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H2A.ub.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H2Aubcomp <- H2Aub[rdTP,]
H2Aub <- H2Aub[onTP,]
M1 <- colMeans(H2Aubcomp)


for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H2Aub_1 <- H2Aub[ind,]
  M1 <- colMeans(H2Aub_1)
  #V <- apply(H2Aub_1, c(2), var) 
  #L <- apply(H2Aub_1, c(2), lowerCI) 
 # U <- apply(H2Aub_1, c(2), upperCI)   
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
 # D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
 # p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H2Aub_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K9me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K9me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K9me3comp <- H3K9me3[rdTP,]
H3K9me3 <- H3K9me3[onTP,]
M1 <- colMeans(H3K9me3comp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K9me3_1 <- H3K9me3[ind,]
  M1 <- colMeans(H3K9me3_1)
  #V <- apply(H3K9me3_1, c(2), var) 
  #L <- apply(H3K9me3_1, c(2), lowerCI) 
  #U <- apply(H3K9me3_1, c(2), upperCI)   
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K9me3_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}


#Now plot an example gene ...
H3K4me1 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K4me1.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me1comp <- H3K4me1[rdTP,]
H3K4me1 <- H3K4me1[onTP,]
M1 <- colMeans(H3K4me1comp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me1_1 <- H3K4me1[ind,]
  M1 <- colMeans(H3K4me1_1)
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  #V <- apply(H3K4me1_1, c(2), var) 
  #L <- apply(H3K4me1_1, c(2), lowerCI) 
  #U <- apply(H3K4me1_1, c(2), upperCI)   
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K4me1_ONTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}


ONList1 <- unique(as.character(AllG[ which(clusters$cluster==1) ,2]))
ONList2 <- unique(as.character(AllG[ which(clusters$cluster==2) ,2]))
ONList3 <- unique(as.character(AllG[ which(clusters$cluster==3) ,2]))
ONList4 <- unique(as.character(AllG[ which(clusters$cluster==4) ,2]))
ONList5 <- unique(as.character(AllG[ which(clusters$cluster==5) ,2]))
ONList6 <- unique(as.character(AllG[ which(clusters$cluster==6) ,2]))
ONList7 <- unique(as.character(AllG[ which(clusters$cluster==7) ,2]))
ONList8 <- unique(as.character(AllG[ which(clusters$cluster==8) ,2]))

saveext <- "./"
write.csv(ONList1, file=paste(saveext,"/Endo_Heatmaps/ONList1.csv",sep=""))
write.csv(ONList2, file=paste(saveext,"/Endo_Heatmaps/ONList2.csv",sep=""))
write.csv(ONList3, file=paste(saveext,"/Endo_Heatmaps/ONList3.csv",sep=""))
write.csv(ONList4, file=paste(saveext,"/Endo_Heatmaps/ONList4.csv",sep=""))
write.csv(ONList5, file=paste(saveext,"/Endo_Heatmaps/ONList5.csv",sep=""))
write.csv(ONList6, file=paste(saveext,"/Endo_Heatmaps/ONList6.csv",sep=""))
write.csv(ONList7, file=paste(saveext,"/Endo_Heatmaps/ONList7.csv",sep=""))
write.csv(ONList8, file=paste(saveext,"/Endo_Heatmaps/ONList8.csv",sep=""))


for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl1.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl2.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl3.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl4.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl5.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl6.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl7.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl8.pdf")

#Cluster by value only
mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]

colnames(mean1) <- labs
NCl <- 4
clusters <- kmeans(mean1, NCl)
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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl_nodelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl1_nodelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl2_nodelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl3_nodelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl4_nodelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl5_nodelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl6_nodelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl7_nodelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl8_nodelta.pdf")





for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl_nodelta_showdelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl1_nodelta_showdelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl2_nodelta_showdelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl3_nodelta_showdelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl4_nodelta_showdelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl5_nodelta_showdelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl6_nodelta_showdelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl7_nodelta_showdelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/ONRDPCl8_nodelta_showdelta.pdf")





#mean2 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONmiusRDTF.pdf")

#mean2 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONminusPermTF.pdf")
##colnames(mean3) <- labs

#Now do off
#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoFull/Up_attributions_OffTP.npy")
mat1 <- np$load("EndoFull/Off_attributions_OffTP.npy")
mat2 <- np$load("EndoFull/Off_attributions_OffTP_perm.npy")

#Index for TP ON memory genes (to get get IDS etc)
offTP <- np$load("EndoFull/OffTP.npy")
offRU <- np$load("EndoFull/OffRU.npy")

#Get the gene IDS etc
Data <-read.table("EndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- as.data.frame(Data$V4[offTP])
colnames(genes) <- c("ID")

AllG <- merge(x = genes, y = TFList, by = "ID", all.x = TRUE)



#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
#mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFF.pdf")

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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl_delta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl1_delta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFNRUPCl2_delta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl3_delta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl4_delta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl5_delta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl6_delta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl7_delta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRDUCl8_delta.pdf")


OFFList1 <- unique(as.character(AllG[ which(clusters$cluster==1) ,2]))
OFFList2 <- unique(as.character(AllG[ which(clusters$cluster==2) ,2]))
OFFList3 <- unique(as.character(AllG[ which(clusters$cluster==3) ,2]))
OFFList4 <- unique(as.character(AllG[ which(clusters$cluster==4) ,2]))
OFFList5 <- unique(as.character(AllG[ which(clusters$cluster==5) ,2]))
OFFList6 <- unique(as.character(AllG[ which(clusters$cluster==6) ,2]))
OFFList7 <- unique(as.character(AllG[ which(clusters$cluster==7) ,2]))
OFFList8 <- unique(as.character(AllG[ which(clusters$cluster==8) ,2]))

saveext <- "./"
write.csv(OFFList1, file=paste(saveext,"/Endo_Heatmaps/OFFList1.csv",sep=""))
write.csv(OFFList2, file=paste(saveext,"/Endo_Heatmaps/OFFList2.csv",sep=""))
write.csv(OFFList3, file=paste(saveext,"/Endo_Heatmaps/OFFList3.csv",sep=""))
write.csv(OFFList4, file=paste(saveext,"/Endo_Heatmaps/OFFList4.csv",sep=""))
write.csv(OFFList5, file=paste(saveext,"/Endo_Heatmaps/OFFList5.csv",sep=""))
write.csv(OFFList6, file=paste(saveext,"/Endo_Heatmaps/OFFList6.csv",sep=""))
write.csv(OFFList7, file=paste(saveext,"/Endo_Heatmaps/OFFList7.csv",sep=""))
write.csv(OFFList8, file=paste(saveext,"/Endo_Heatmaps/OFFList8.csv",sep=""))

#Now plot an example gene ...
H3K4me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K4me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me3comp <- H3K4me3[ruTP,]
H3K4me3 <- H3K4me3[offTP,]
M1 <- colMeans(H3K4me3comp)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me3_1 <- H3K4me3[ind,]
  M <- colMeans(H3K4me3_1)
  #V <- apply(H3K4me3_1, c(2), var) 
  #L <- apply(H3K4me3_1, c(2), lowerCI) 
  #U <- apply(H3K4me3_1, c(2), upperCI)   
  
  D <- data.frame(X=c( seq(1,length(M),1), seq(1,length(M1),1))
                  ,M=c(M,M1),type=as.factor(c(rep(0, length(M)), rep(1, length(M1))) ) )
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M,color = type), size = 1)+ theme_bw()
  
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K4me3_OFfTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K36me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K36me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K36me3comp <- H3K36me3[ruTP,]
H3K36me3 <- H3K36me3[offTP,]
M1 <- colMeans(H3K36me3_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K36me3_1 <- H3K36me3[ind,]
  M <- colMeans(H3K36me3_1)
  #V <- apply(H3K36me3_1, c(2), var) 
  #L <- apply(H3K36me3_1, c(2), lowerCI) 
  #U <- apply(H3K36me3_1, c(2), upperCI)   
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K36me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K27ac <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K27ac.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K27accomp <- H3K27ac[ruTP,]
H3K27ac <- H3K27ac[offTP,]
M1 <- colMeans(H3K27ac_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K27ac_1 <- H3K27ac[ind,]
  M <- colMeans(H3K27ac_1)
  #V <- apply(H3K27ac_1, c(2), var) 
  #L <- apply(H3K27ac_1, c(2), lowerCI) 
  #U <- apply(H3K27ac_1, c(2), upperCI)   
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K27ac_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K79me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K79me3.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K79me3comp <- H3K79me3[ruTP,]
H3K79me3 <- H3K79me3[offTP,]
M1 <- colMeans(H3K79me3_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K79me3_1 <- H3K79me3[ind,]
  M <- colMeans(H3K79me3_1)
  #V <- apply(H3K79me3_1, c(2), var) 
  #L <- apply(H3K79me3_1, c(2), lowerCI) 
  #U <- apply(H3K79me3_1, c(2), upperCI)   
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K79me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H2Aub <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H2A.ub.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H2Aubcomp <- H2Aub[ruTP,]
H2Aub <- H2Aub[offTP,]
M1 <- colMeans(H2Aub_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H2Aub_1 <- H2Aub[ind,]
  M <- colMeans(H2Aub_1)
  #V <- apply(H2Aub_1, c(2), var) 
  #L <- apply(H2Aub_1, c(2), lowerCI) 
  #U <- apply(H2Aub_1, c(2), upperCI)   
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H2Aub_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K9me3 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K9me3.merged.sorted.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K9me3 <- H3K9me3[offTP,]
H3K9me3comp <- H3K9me3[ruTP,]
M1 <- colMeans(H3K9me3_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K9me3_1 <- H3K9me3[ind,]
  M <- colMeans(H3K9me3_1)
  #V <- apply(H3K9me3_1, c(2), var) 
  #L <- apply(H3K9me3_1, c(2), lowerCI) 
  #U <- apply(H3K9me3_1, c(2), upperCI)   
  
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  #p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K9me3_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}

#Now plot an example gene ...
H3K4me1 <- read.table("/Volumes/Overflow/HistoneMatrix/Endo/Endo_H3K4me1.singleRep.bam_vs_Input_CPM.bwe.tab", skip = 3)
H3K4me1comp <- H3K4me1[ruTP,]
H3K4me1 <- H3K4me1[offTP,]
M1 <- colMeans(H3K4me1_1)

for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  H3K4me1_1 <- H3K4me1[ind,]
  M <- colMeans(H3K4me1_1)
  #V <- apply(H3K4me1_1, c(2), var) 
  #L <- apply(H3K4me1_1, c(2), lowerCI) 
  #U <- apply(H3K4me1_1, c(2), upperCI)   
  #D <- data.frame(X=seq(1,length(M),1),M=M,U=U,L=L)
  p <- ggplot(D, aes(X, M))+geom_line(aes(x=X, y=M), size = 1)+geom_ribbon(aes(x=X, ymin=L, ymax=U), linetype=2, alpha=0.1)+ theme_bw()
  ggsave(filename=paste("Endo_Heatmaps/H3K4me1_OFFTPCl",i,".pdf",sep=""),width = 16, height = 8, p)
}





for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl1.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFNRUPCl2.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl3.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl4.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl5.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl6.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl7.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRDUCl8.pdf")










mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean2 = apply(mat1-mat, c(2,3), sum) / dim(mat1)[1]

colnames(mean1) <- labs
NCl <- 4
clusters <- kmeans(mean1, NCl)
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
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl_nodelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl1_nodelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFNRUPCl2_nodelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl3_nodelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl4_nodelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl5_nodelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl6_nodelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl7_nodelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRDUCl8_nodelta.pdf")


for (i in 1:NCl)  {
  ind <- which(clusters$cluster==i)
  uniqueMat1[i,,] <- apply(mat1[ind,,]-mat[ind,,], c(2,3), sum) / length(ind)
  SummaryMat1[i,] <- 600*colSums(mean1[ind,])/length(ind)
}
colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- c(paste("Cl",rep(1:NCl, each=1)) )

mat_breaks <- seq(-0.2, 0.2, length.out = 20)
pheatmap(t(SummaryMat1), breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl_nodelta_showdelta.pdf")

mat_breaks <- seq(-0.2/100, 0.2/100, length.out = 20)

x <- t(uniqueMat1[1,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl1_nodelta_showdelta.pdf")

x <- t(uniqueMat1[2,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFNRUPCl2_nodelta_showdelta.pdf")

x <- t(uniqueMat1[3,,])
rownames(x) <- labs
pheatmap( x,   breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl3_nodelta_showdelta.pdf")

x <- t(uniqueMat1[4,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl4_nodelta_showdelta.pdf")

x <- t(uniqueMat1[5,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl5_nodelta_showdelta.pdf")

x <- t(uniqueMat1[6,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl6_nodelta_showdelta.pdf")

x <- t(uniqueMat1[7,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRUPCl7_nodelta_showdelta.pdf")

x <- t(uniqueMat1[8,,])
rownames(x) <- labs
pheatmap( x,  breaks = mat_breaks, color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "Endo_Heatmaps/OFFRDUCl8_nodelta_showdelta.pdf")


#Order of experiments
#labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
#redblue1 <- colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))


#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
#mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
#mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFF.pdf")

#mean2 = apply(mat1[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat1[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat1[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "OFFTF.pdf")

#mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
#mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RU.pdf")

#mean2 = apply(mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RUTF.pdf")


#mean2 = apply(mat1-mat, c(1,3), sum) / dim(mat)[2]
#mean3 = apply(mat1-mat, c(2,3), sum) / dim(mat)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFmiusRU.pdf")

#mean2 = apply(mat1-mat2, c(1,3), sum) / dim(mat)[2]
#mean3 = apply(mat1-mat2, c(2,3), sum) / dim(mat)[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFminusPerm.pdf")

#mean2 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,]-mat[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFmiusRUTF.pdf")

#mean2 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(1,3), sum) / dim(mat[which(AllG[,5]==1),,])[2]
#mean3 = apply(mat1[which(AllG[,5]==1),,]-mat2[which(AllG[,5]==1),,], c(2,3), sum) / dim(mat[which(AllG[,5]==1),,])[1]
#colnames(mean3) <- labs
#pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "OFFminusPermTF.pdf")


