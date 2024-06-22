#Innstall reticulate (uncommet the line below to do so)
#install.packages("reticulate")
library(reticulate)
library(pheatmap)

setwd("/Volumes/Overflow")
np <- import("numpy")

#Load activtion data in this case the RD activations for the TP ON memory genes
mat <- np$load("EndoFull/Down_attributions_OnTP.npy")
mat1 <- np$load("EndoFull/On_attributions_OnTP.npy")
mat2 <- np$load("EndoFull/On_attributions_OnTP_perm.npy")

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
mean2 = apply(mat, c(1,3), sum) / dim(mat)[2]
mean3 = apply(mat, c(2,3), sum) / dim(mat)[1]
colnames(mean3) <- labs
pheatmap(t(mean3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, main = c("Average over all genes"))

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
NCl <- 15
clusters <- kmeans(mean2, NCl)
SummaryMat <- array(0,c(NCl,dim(mat)[3]))

for (i in 1:NCl)  {
    ind <- which(clusters$cluster==i)
    uniqueMat[i,,] <- apply(mat[ind,,], c(2,3), sum) / length(idl)
    #genes_cl <- genes[ind] pull out the genes that belong, can write to csv if need be
    #pheatmap(t(M3),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Subset of genes"))
    SummaryMat[i,] <- colSums(mean2[ind,])/length(ind)
}

colnames(SummaryMat) <- labs
rownames(SummaryMat) <- paste("Cl",rep(1:NCl, each=1))
pheatmap(t(SummaryMat),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=TRUE, main = c("Cluster summary"))

#Look at specific cluster
cltolook <- 5
Msp <- uniqueMat[cltolook,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"))s

#Look at specific cluster
cltolook <- 3
Msp <- uniqueMat[cltolook,,]
colnames(Msp) <- labs
pheatmap(t(Msp),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, main = c("Cluster"))

##### What about the delta's?
#Can load in a seperate set of activatios to calculate the delta.
mat2 <- np$load("EndoFull/NUMPY/ON_attributions_OnTP.npy")
Delta_mat = mat2 - mat #On - Down @ OnTP

pheatmap(t(Delta_mat[genechoice,,]),color =  redblue1(20),  gaps_row=c(10,20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE)


#Now should roll this out. Need to look at the activations, deltas, and clusterings
library(ggplot2)
TPR <- np$load("EndoFull/NUMPY/tpr_kearst.npy")
FPR <- np$load("EndoFull/NUMPY/fpr_kearst.npy")

AUC <- data.frame(TPR,FPR)

ggplot(data=AUC, aes(x=fpr, y=tpr)) + geom_line() + theme_minimal() + xlim(0,1)+ylim(0,1)
#ggsave(filename=paste("~/Desktop/NoGenesByStage.pdf",sep=""), plot = p, width = 16, height = 16)
