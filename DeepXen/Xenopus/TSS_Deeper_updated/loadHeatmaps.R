library(reticulate)


np <- import("numpy")

#Load data

mat <- np$load("/Users/christopherpenfold/Desktop/Code/DeepXen/DeepXen/Xenopus/TSS_Deeper_updated/EndoFullSomEndo/hm1_delta.npy")

labs<-c("Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation")


redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(mat,color =  redblue1(20),  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/HM5.pdf")
