#Innstall reticulate (uncommet the line below to do so)
#install.packages("reticulate")
library(reticulate)
library(pheatmap)

setwd("/Volumes/Overflow")
np <- import("numpy")

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoChr/Down_attributions_OnTP.npy")
mat1 <- np$load("EndoChr/On_attributions_OnTP.npy")
mat2 <- np$load("EndoChr/On_attributions_OnTP_perm.npy")

#mat1 <- np$load("EndoFull/Down_attributions_OnTP_1.npy")


#Index for TP ON memory genes (to get get IDS etc)
onTP <- np$load("EndoFull/OnTP.npy")

#Get the gene IDS etc
Data <-read.table("EndoFull/allpeaks_labelled_cum_final.se.bed",sep="\t", header = F, row.names=NULL)
genes <- Data$V4[onTP]

#Order of experiments
labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")

#Colours for the heatmap
redblue1 <- colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))

#Example to print out the heatmap and save it
#mat_breaks <- seq(0, 2, length.out = 20) #Set the scale bar manually rather than automatically
#pheatmap(mat,color =  redblue1(20),  border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE,  filename = "HMexample.pdf")

#1) First plot an example heatmap e.g., just pick one sample and show histone mod by location
genechoice <- 1
M <- t(mat[genechoice,,])
rownames(M) <- labs
mat_breaks <- seq(0, 0.01, length.out = 20) #Set the scale bar manually rather than automatically
genechoice <- 22
pheatmap(M,color =  redblue1(20),  gaps_row=c(10,20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = as.character(genes[genechoice]))

#M2 <- t(mat1[genechoice,,])
#colnames(M2) <- c('Donor','Target')
#pheatmap(M2,color =  redblue1(20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE)

#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean2 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean3 = apply(mat1, c(2,3), sum) / dim(mat1)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "ON.pdf")

mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"), filename = "RD.pdf")

mean2 = apply(mat1-mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONmiusRD.pdf")

mean2 = apply(mat1-mat2, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat1-mat2, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"),  filename = "ONminusPerm.pdf")


#3) Can prefilter this by a list. In this case we will just generate a random list.
rownames(mean2) <- genes
colnames(mean2) <- labs
slist <- genes[seq(500, 800, 10)]

M3 <- mean2[slist,]
pheatmap(t(M3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Subset of genes"))

#Try looking at just a single mark for all genes in the list. Note this will take a while.
inde  <- which(labs=="Endo_H3K4me3")
M4 <- mat[,,inde]
#pheatmap(t(M4),color =  redblue1(20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Subset of marks"),filename = "HMexample.pdf")

#Unique genes only (average over the others)
#uniqueMat <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
#for (i in 1:length(unique(genes)) ) {
#    idl <- which(genes==unique(genes)[i])
#
#    if (length(idl)==1) {
#    uniqueMat[i,,] <- mat[idl,,]
#    } else {
#    uniqueMat[i,,] <- apply(mat[idl,,], c(2,3), sum) / length(idl)
#    }
#}
#twoDimMat <- matrix(uniqueMat, dim(uniqueMat)[1], dim(uniqueMat)[2]*dim(uniqueMat)[3])
#rownames(twoDimMat) <- unique(genes)

#Do k-means on average histone values and then look  at the overall profiles.
#2) Now average over all genes in this set. This set can be a seperate list loaded in of e.g., histone modifiers etc.
mean1 = apply(mat1, c(1,3), sum) / dim(mat1)[2]
mean2 = apply(mat1, c(2,3), sum) / dim(mat1)[1]

mean3 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean4 = apply(mat, c(2,3), sum) / dim(mat)[1]

mean5 = apply(mat2, c(1,3), sum) / dim(mat2)[2]
mean6 = apply(mat2, c(2,3), sum) / dim(mat2)[1]

colnames(mean3) <- labs
NCl <- 15
clusters <- kmeans(mean1, NCl)
SummaryMat1 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat2 <- array(0,c(NCl,dim(mat)[3]))
SummaryMat3 <- array(0,c(NCl,dim(mat)[3]))
uniqueMat1 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat2 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))
uniqueMat3 <- array(0,c(length(unique(genes)),dim(mat)[2],dim(mat)[3]))

for (i in 1:NCl)  {
    ind <- which(clusters$cluster==i)
    uniqueMat1[i,,] <- apply(mat1[ind,,], c(2,3), sum) / length(ind)
    uniqueMat2[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(ind)
    uniqueMat3[i,,] <- apply(mat2[ind,,], c(2,3), sum) / length(ind)
    #genes_cl <- genes[ind] pull out the genes that belong, can write to csv if need be
    #pheatmap(t(M3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Subset of genes"))
    SummaryMat1[i,] <- colSums(mean1[ind,])/length(ind)
    SummaryMat2[i,] <- colSums(mean3[ind,])/length(ind)
    SummaryMat3[i,] <- colSums(mean5[ind,])/length(ind)
}
SummaryMat <- rbind(SummaryMat1,SummaryMat2,SummaryMat3)
colnames(SummaryMat) <- labs
#colnames(SummaryMat) <- labs
rownames(SummaryMat) <- c(paste("Cl",rep(1:NCl, each=1)), paste("ClR",rep(1:NCl, each=1)), paste("ClP",rep(1:NCl, each=1)) )
pheatmap(t(SummaryMat),color =  redblue1(20), gaps_col = c(NCl,NCl*2),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster summary"),filename = "ONRDPCl.pdf")

colnames(SummaryMat1) <- labs
rownames(SummaryMat1) <- paste("Cl",rep(1:NCl, each=1))
pheatmap(t(SummaryMat1),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"),filename = "ONCl.pdf")


#Look at specific cluster
for (i in seq(1,NCl,by=1)) {
Msp <- uniqueMat1[i,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"),filename = paste("ONCl", i, ".pdf", sep="") )

Msp <- uniqueMat2[i,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"),filename = paste("RDCl", i, ".pdf", sep="") )

Msp <- uniqueMat1[i,,]-uniqueMat2[i,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"),filename = paste("ONClminusRD", i, ".pdf", sep="") )

Msp <- uniqueMat1[i,,]-uniqueMat3[i,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"),filename = paste("ONClminusPerm", i, ".pdf", sep="") )

}
#Look at specific cluster
cltolook <- 3
Msp <- uniqueMat1[cltolook,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"))

##### What about the delta's?
#Can load in a seperate set of activatios to calculate the delta.
#mat2 <- np$load("EndoFull/NUMPY/ON_attributions_OnTP.npy")
#Delta_mat = mat2 - mat #On - Down @ OnTP

#pheatmap(t(Delta_mat[genechoice,,]),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE)


#Now should roll this out. Need to look at the activations, deltas, and clusterings
#library(ggplot2)
#TPR <- np$load("EndoFull/NUMPY/tpr_kearst.npy")
#FPR <- np$load("EndoFull/NUMPY/fpr_kearst.npy")

#AUC <- data.frame(TPR,FPR)

#ggplot(data=AUC, aes(x=fpr, y=tpr)) + geom_line() + theme_minimal() + xlim(0,1)+ylim(0,1)
#ggsave(filename=paste("~/Desktop/NoGenesByStage.pdf",sep=""), plot = p, width = 16, height = 16)
