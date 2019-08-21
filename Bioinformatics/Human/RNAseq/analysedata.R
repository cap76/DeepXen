
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")
library(dplyr)
library(ggplot2)
library("ggbeeswarm")
library("RColorBrewer")
library("gplots")
library(GenomicFeatures)


#sampleNames <- c("HDF_1","HDF_2","HDF_3","IVF1_CM1","IVF1_CM2","IVF1_EC1","IVF1_EC2","IVF1_ESC1","IVF1_ESC2","IVF2_CM1","IVF2_CM2","IVF2_EC1","IVF2_EC2","IVF2_ESC1","IVF2_ESC2","NT1_CM1","NT1_CM2","NT1_EC1","NT1_EC2","NT1_NT1","NT1_NT2","NT2_CM1","NT2_CM2","NT2_EC1","NT2_EC2","NT2_NT1","NT2_NT2")
#
#"IVF1_ESC","IVF1_ESC","IVF2_ESC","IVF2_ESC","NT1_ESC","NT1_ESC","NT2_ESC","NT2_ESC"
#SRR5217330=ESO7
#331=ESo7CM
#332 ES07 EC
#33 EC
#34 ESC
#35 ESC
#SRR5217336=ESO8
#337=ESo7CM
#338 ES07 EC
#339 EC
#340 ESC
#341 ESC
#42  i12C-CM
#43  CM
#44  EC
#45  EC
#46  i12C-CM
#47
#48 i12J-CM_rep1
#49
#50
#51



sampleNames <- c("HDF_1","HDF_2","HDF_3","IVF1_CM1","IVF1_CM2","IVF1_EC1","IVF1_EC2","IVF1_ESC1","IVF1_ESC2","IVF2_CM1","IVF2_CM2","IVF2_EC1","IVF2_EC2","IVF2_ESC1","IVF2_ESC2","NT1_CM1","NT1_CM2","NT1_EC1","NT1_EC2","NT1_NT1","NT1_NT2","NT2_CM1","NT2_CM2","NT2_EC1","NT2_EC2","NT2_NT1","NT2_NT2","iPSC1_CM1","iPSC1_CM2","iPSC1_EC1","iPSC1_EC2","iPSC1_PSC1","iPSC1_PSC2","iPSC2_CM1","iPSC2_CM2","iPSC2_EC1","iPSC2_EC2","iPSC2_PSC1","iPSC2_PSC2")
sampleFiles <- c("SRR5029372/SRR5029372.counts","SRR5029373/SRR5029373.counts","SRR5029374/SRR5029374.counts","SRR5217330/SRR5217330.counts","SRR5217331/SRR5217331.counts","SRR5217332/SRR5217332.counts","SRR5217333/SRR5217333.counts","SRR5217334/SRR5217334.counts","SRR5217335/SRR5217335.counts","SRR5217336/SRR5217336.counts","SRR5217337/SRR5217337.counts","SRR5217338/SRR5217338.counts","SRR5217339/SRR5217339.counts","SRR5217340/SRR5217340.counts","SRR5217341/SRR5217341.counts","SRR5217354/SRR5217354.counts","SRR5217355/SRR5217355.counts","SRR5217356/SRR5217356.counts","SRR5217357/SRR5217357.counts","SRR5217358/SRR5217358.counts","SRR5217359/SRR5217359.counts","SRR5217360/SRR5217360.counts","SRR5217361/SRR5217361.counts","SRR5217362/SRR5217362.counts","SRR5217363/SRR5217363.counts","SRR5217364/SRR5217364.counts","SRR5217365/SRR5217365.counts","SRR5217342/SRR5217342.counts","SRR5217343/SRR5217343.counts","SRR5217344/SRR5217344.counts","SRR5217345/SRR5217345.counts","SRR5217346/SRR5217346.counts","SRR5217347/SRR5217347.counts","SRR5217348/SRR5217348.counts","SRR5217349/SRR5217349.counts","SRR5217350/SRR5217350.counts","SRR5217351/SRR5217351.counts","SRR5217352/SRR5217352.counts","SRR5217353/SRR5217353.counts")
sampleCondition <- c("HDF","HDF","HDF","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT","iPSC_CM","iPSC_CM","iPSC_EC","iPSC_EC","iPSC_PSC","iPSC_PSC","iPSC_CM","iPSC_CM","iPSC_EC","iPSC_EC","iPSC_PSC","iPSC_PSC")
sampleReplicate <- c("1","2","3","1","2","1","2","1","2","1","2","1","2","1","2","3","4","3","4","3","4","3","4","3","4","3","4","1","2","1","2","1","2","3","4","3","4","3","4")

txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
geneLength <- as.data.frame(unlist(exonic.gene.sizes))




#"SRR2501664/SRR2501664.counts"



#sampleFiles <- c("SRR5029372/SRR5029372.counts","SRR5029373/SRR5029373.counts","SRR5029374/SRR5029374.counts","SRR5217330/SRR5217330.counts","SRR5217331/SRR5217331.counts","SRR5217332/SRR5217332.counts","SRR5217333/SRR5217333.counts","SRR5217334/SRR5217334.counts","SRR5217335/SRR5217335.counts","SRR5217336/SRR5217336.counts","SRR5217337/SRR5217337.counts","SRR5217338/SRR5217338.counts","SRR5217339/SRR5217339.counts","SRR5217340/SRR5217340.counts","SRR5217341/SRR5217341.counts","SRR5217354/SRR5217354.counts","SRR5217355/SRR5217355.counts","SRR5217356/SRR5217356.counts","SRR5217357/SRR5217357.counts","SRR5217358/SRR5217358.counts","SRR5217359/SRR5217359.counts","SRR5217360/SRR5217360.counts","SRR5217361/SRR5217361.counts","SRR5217362/SRR5217362.counts","SRR5217363/SRR5217363.counts","SRR5217364/SRR5217364.counts","SRR5217365/SRR5217365.counts")
#sampleCondition <- c("HDF","HDF","HDF","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT")
#sampleReplicate <- c("1","2","3","1","2","1","2","1","2","3","4","3","4","3","4","1","2","1","2","1","2","3","4","3","4","3","4")

# Set the working directory
directory <- "."
setwd(directory)

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, rep = sampleReplicate)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/mnt/scratch/gurdon/cap76/DeepXen/Human/RNAseq/", design=~condition)

#Now we get gene lengths
txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
geneLength <- as.data.frame(unlist(exonic.gene.sizes))

#Extract out useful information
dds <- DESeq(ddsHTSeq)
res <- results(dds)
fpms <- fpm(dds)

#Now generate the 
y <- merge(x = fpms, y = geneLength, by=0, all.x = TRUE)
colnames(y)[length(colnames(y))] <- "basepairs"

mcols(dds) <- cbind(mcols(dds), y)

fpkms <- fpkm(dds)

#HDF vs others
res1<-results(dds, contrast=c("condition","HDF","IVF1_CM"))
res2<-results(dds, contrast=c("condition","HDF","IVF1_EC"))
res3<-results(dds, contrast=c("condition","HDF","IVF2_CM"))
res4<-results(dds, contrast=c("condition","HDF","IVF2_EC"))

res5<-results(dds, contrast=c("condition","HDF","NT1_CM"))
res6<-results(dds, contrast=c("condition","HDF","NT1_EC"))
res7<-results(dds, contrast=c("condition","HDF","NT2_CM"))
res8<-results(dds, contrast=c("condition","HDF","NT2_EC"))

#CM comparisons
res9<-results(dds, contrast=c("condition","NT1_CM","IVF1_CM"))
res10<-results(dds, contrast=c("condition","NT2_CM","IVF1_CM"))
res11<-results(dds, contrast=c("condition","NT1_CM","IVF2_CM"))
res12<-results(dds, contrast=c("condition","NT2_CM","IVF2_CM"))

#EC comparisons
res13<-results(dds, contrast=c("condition","NT1_EC","IVF1_EC"))
res14<-results(dds, contrast=c("condition","NT1_EC","IVF2_EC"))
res15<-results(dds, contrast=c("condition","NT2_EC","IVF1_EC"))
res16<-results(dds, contrast=c("condition","NT2_EC","IVF2_EC"))

#EC comparisons
res13<-results(dds, contrast=c("condition","NT1_EC","IVF1_EC"))
res14<-results(dds, contrast=c("condition","NT1_EC","IVF2_EC"))
res15<-results(dds, contrast=c("condition","NT2_EC","IVF1_EC"))
res16<-results(dds, contrast=c("condition","NT2_EC","IVF2_EC"))

#ESC comparisons
res17<-results(dds, contrast=c("condition","NT1_NT","IVF1_ESC"))
res18<-results(dds, contrast=c("condition","NT2_NT","IVF1_ESC"))
res19<-results(dds, contrast=c("condition","NT1_NT","IVF2_ESC"))
res20<-results(dds, contrast=c("condition","NT2_NT","IVF2_ESC"))

#HDF vs others
res21<-results(dds, contrast=c("condition","HDF","IVF1_ESC"))
res22<-results(dds, contrast=c("condition","HDF","IVF2_ESC"))
res23<-results(dds, contrast=c("condition","HDF","NT1_NT"))
res24<-results(dds, contrast=c("condition","HDF","NT2_NT"))

library(ggplot2)
library(beeswarm)
library(ggbeeswarm)

D1 <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res9$log2FoldChange,res9$padj))
p <- ggplot(D1, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red") 
ggsave("NT1_CM_vs_IVF1.pdf")
dev.off()

D2 <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res10$log2FoldChange,res10$padj))
p <- ggplot(D2, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_CM_vs_IVF1.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res11$log2FoldChange,res11$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT1_CM_vs_IVF2.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res12$log2FoldChange,res12$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_CM_vs_IVF2.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res13$log2FoldChange,res13$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("EC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT1_EC_vs_IVF1.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res14$log2FoldChange,res14$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("EC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_EC_vs_IVF1.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res15$log2FoldChange,res15$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("EC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT1_EC_vs_IVF2.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res16$log2FoldChange,res16$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("EC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_EC_vs_IVF2.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res17$log2FoldChange,res17$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("ESC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT1_ESC_vs_IVF1.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res18$log2FoldChange,res18$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("ESC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_ESC_vs_IVF1.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res19$log2FoldChange,res19$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("ESC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT1_ESC_vs_IVF2.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res20$log2FoldChange,res20$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("ESC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT2_ESC_vs_IVF2.pdf")
dev.off()






#pdf('NT1_CM_vs_IVF1.pdf')
#p
#dev.close()


#rd1<-as.data.frame(res1)
#rd2<-as.data.frame(res2)
#rd3<-as.data.frame(res3)
#rd4<-as.data.frame(res4)
#rd5<-as.data.frame(res5)
#rd6<-as.data.frame(res6)
#rd7<-as.data.frame(res7)
#rd8<-as.data.frame(res8)
#rd9<-as.data.frame(res9)
#rd10<-as.data.frame(res10)
#rd10<-as.data.frame(res11)
#rd10<-as.data.frame(res12)
#rd10<-as.data.frame(res13)
#rd10<-as.data.frame(res14)
#rd10<-as.data.frame(res15)
#rd10<-as.data.frame(res16)
#
#write.csv(as.data.frame(fpkms), file = "FPKM.csv")
#write.csv(as.data.frame(fpms), file = "FPM.csv")
#
ON_HDF_IVF1_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF1_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res18)$pval<0.05 & as.data.frame(res18)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res19)$pval<0.05 & as.data.frame(res19)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res20)$pval<0.05 & as.data.frame(res20)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

ON_HDF_IVF1_NT1_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF1_NT2_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT1_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT2_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

ON_HDF_IVF1_NT1_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF1_NT2_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT1_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res15)$pval<0.05 & as.data.frame(res15)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF2_NT2_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res16)$pval<0.05 & as.data.frame(res16)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

#   #donor vs ivf <0.05 & FC>0, nt vs ivf < 0.05 & FC>0, & fpkm>1
OFF_HDF_IVF1_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange<0)
OFF_HDF_IVF1_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res18)$pval<0.05 & as.data.frame(res18)$log2FoldChange<0)
OFF_HDF_IVF2_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res19)$pval<0.05 & as.data.frame(res19)$log2FoldChange<0)
OFF_HDF_IVF2_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res20)$pval<0.05 & as.data.frame(res20)$log2FoldChange<0)

OFF_HDF_IVF1_NT1_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange<0)
OFF_HDF_IVF1_NT2_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange<0)
OFF_HDF_IVF2_NT1_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange<0)
OFF_HDF_IVF2_NT2_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange<0)

OFF_HDF_IVF1_NT1_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange<0)
OFF_HDF_IVF1_NT2_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange<0)
OFF_HDF_IVF2_NT1_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res15)$pval<0.05 & as.data.frame(res15)$log2FoldChange<0)
OFF_HDF_IVF2_NT2_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res16)$pval<0.05 & as.data.frame(res16)$log2FoldChange<0)

#
RD_HDF_IVF1_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res23)$pval<0.05 & as.data.frame(res23)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF1_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res24)$pval<0.05 & as.data.frame(res24)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res23)$pval<0.05 & as.data.frame(res23)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res24)$pval<0.05 & as.data.frame(res24)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

RD_HDF_IVF1_NT1_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF1_NT2_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT1_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT2_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

RD_HDF_IVF1_NT1_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF1_NT2_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT1_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF2_NT2_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

#   #donor vs ivf <0.05 & FC>0, nt vs ivf < 0.05 & FC>0, & fpkm>1
RU_HDF_IVF1_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res23)$pval<0.05 & as.data.frame(res23)$log2FoldChange<0 & as.data.frame(res17)$pval>0.05)
RU_HDF_IVF1_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res24)$pval<0.05 & as.data.frame(res24)$log2FoldChange<0 & as.data.frame(res18)$pval>0.05)
RU_HDF_IVF2_NT1_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res23)$pval<0.05 & as.data.frame(res23)$log2FoldChange<0 & as.data.frame(res19)$pval>0.05)
RU_HDF_IVF2_NT2_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res24)$pval<0.05 & as.data.frame(res24)$log2FoldChange<0 & as.data.frame(res20)$pval>0.05)

RU_HDF_IVF1_NT1_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange<0 & as.data.frame(res9)$pval>0.05)
RU_HDF_IVF1_NT2_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange<0 & as.data.frame(res10)$pval>0.05)
RU_HDF_IVF2_NT1_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange<0 & as.data.frame(res11)$pval>0.05)
RU_HDF_IVF2_NT2_CM <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange<0 & as.data.frame(res12)$pval>0.05)

RU_HDF_IVF1_NT1_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange<0 & as.data.frame(res13)$pval>0.05)
RU_HDF_IVF1_NT2_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange<0 & as.data.frame(res14)$pval>0.05)
RU_HDF_IVF2_NT1_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange<0 & as.data.frame(res15)$pval>0.05)
RU_HDF_IVF2_NT2_EC <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange<0 & as.data.frame(res16)$pval>0.05)





#Merge the cell lines!
sampleNames <- c("HDF_1","HDF_2","HDF_3","IVF1_CM1","IVF1_CM2","IVF1_EC1","IVF1_EC2","IVF1_ESC1","IVF1_ESC2","IVF2_CM1","IVF2_CM2","IVF2_EC1","IVF2_EC2","IVF2_ESC1","IVF2_ESC2","NT1_CM1","NT1_CM2","NT1_EC1","NT1_EC2","NT1_NT1","NT1_NT2","NT2_CM1","NT2_CM2","NT2_EC1","NT2_EC2","NT2_NT1","NT2_NT2")
sampleFiles <- c("SRR5029372/SRR5029372.counts","SRR5029373/SRR5029373.counts","SRR5029374/SRR5029374.counts","SRR5217330/SRR5217330.counts","SRR5217331/SRR5217331.counts","SRR5217332/SRR5217332.counts","SRR5217333/SRR5217333.counts","SRR5217334/SRR5217334.counts","SRR5217335/SRR5217335.counts","SRR5217336/SRR5217336.counts","SRR5217337/SRR5217337.counts","SRR5217338/SRR5217338.counts","SRR5217339/SRR5217339.counts","SRR5217340/SRR5217340.counts","SRR5217341/SRR5217341.counts","SRR5217354/SRR5217354.counts","SRR5217355/SRR5217355.counts","SRR5217356/SRR5217356.counts","SRR5217357/SRR5217357.counts","SRR5217358/SRR5217358.counts","SRR5217359/SRR5217359.counts","SRR5217360/SRR5217360.counts","SRR5217361/SRR5217361.counts","SRR5217362/SRR5217362.counts","SRR5217363/SRR5217363.counts","SRR5217364/SRR5217364.counts","SRR5217365/SRR5217365.counts")
sampleCondition <- c("HDF","HDF","HDF","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","IVF_CM","IVF_CM","IVF_EC","IVF_EC","IVF_ESC","IVF_ESC","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT","NT_CM","NT_CM","NT_EC","NT_EC","NT_NT","NT_NT")
sampleReplicate <- c("1","2","3","1","2","1","2","1","2","1","2","1","2","1","2","1","2","1","2","1","2","1","2","1","2","1","2")

# Set the working directory
directory <- "."
setwd(directory)

sampleTable2 <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, rep = sampleReplicate)

ddsHTSeq2<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable2, directory="/mnt/scratch/gurdon/cap76/DeepXen/Human/RNAseq/", design=~condition)

#Now we get gene lengths
txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
geneLength <- as.data.frame(unlist(exonic.gene.sizes))

#Extract out useful information
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
fpms2 <- fpm(dds2)

#Now generate the
y2 <- merge(x = fpms2, y = geneLength, by=0, all.x = TRUE)
colnames(y2)[29] <- "basepairs"

mcols(dds2) <- cbind(mcols(dds2), y2)

fpkms <- fpkm(dds2)

#HDF vs others
res1<-results(dds, contrast=c("condition","HDF","IVF_CM"))
res21<-results(dds2, contrast=c("condition","HDF","IVF_ESC"))
res2<-results(dds2, contrast=c("condition","HDF","IVF_EC"))
res5<-results(dds2, contrast=c("condition","HDF","NT_CM"))
res6<-results(dds2, contrast=c("condition","HDF","NT_EC"))
res10<-results(dds2, contrast=c("condition","HDF","NT_NT"))
res7<-results(dds2, contrast=c("condition","HDF","iPSC_CM"))
res8<-results(dds2, contrast=c("condition","HDF","iPSC_EC"))
res9<-results(dds2, contrast=c("condition","HDF","iPSC_PSC"))
#NT comparisons
res9<-results(dds2, contrast=c("condition","NT_CM","IVF_CM"))
res13<-results(dds2, contrast=c("condition","NT_EC","IVF_EC"))
res17<-results(dds2, contrast=c("condition","NT_NT","IVF_ESC"))
#iPSC
res11<-results(dds2, contrast=c("condition","iPSC_CM","IVF_CM"))
res12<-results(dds2, contrast=c("condition","iPSC_EC","IVF_EC"))
#ESC comparisons
res14<-results(dds2, contrast=c("condition","iPSC_PSC","IVF_ESC"))


ON_HDF_IVF_NT_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF_NT_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF_NT_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

ON_HDF_IVF_iPSC_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF_iPSC_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
ON_HDF_IVF_iPSC_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

OFF_HDF_IVF_NT_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange<0)
OFF_HDF_IVF_NT_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange<0)
OFF_HDF_IVF_NT_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange<0)

OFF_HDF_IVF_iPSC_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange<0)
OFF_HDF_IVF_iPSC_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange<0)
OFF_HDF_IVF_iPSC_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange<0)


RD_HDF_IVF_NT_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF_NT_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF_NT_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

RD_HDF_IVF_iPSC_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange>0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF_iPSC_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RD_HDF_IVF_iPSC_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)

RU_HDF_IVF_NT_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange<0 & as.data.frame(res17)$pval>0.05)
RU_HDF_IVF_NT_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange<0 & as.data.frame(res9)$pval>0.05)
RU_HDF_IVF_NT_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange<0 & as.data.frame(res13)$pval>0.05)

RU_HDF_IVF_iPSC_ESC <- which(as.data.frame(res21)$pval<0.05 & as.data.frame(res21)$log2FoldChange<0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange<0 & as.data.frame(res14)$pval>0.05)
RU_HDF_IVF_iPSC_CM <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange<0 & as.data.frame(res11)$pval>0.05)
RU_HDF_IVF_iPSC_EC <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange<0 & as.data.frame(res12)$pval>0.05)



write.table(rownames(res21)[ON_HDF_IVF_NT_ESC], file = "ON_HDF_IVF_NT_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[ON_HDF_IVF_NT_CM], file = "ON_HDF_IVF_NT_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[ON_HDF_IVF_NT_EC], file = "ON_HDF_IVF_NT_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[ON_HDF_IVF_iPSC_ESC], file = "ON_HDF_IVF_iPSC_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[ON_HDF_IVF_iPSC_CM], file = "ON_HDF_IVF_iPSC_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[ON_HDF_IVF_iPSC_EC], file = "ON_HDF_IVF_iPSC_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[OFF_HDF_IVF_NT_ESC], file = "OFF_HDF_IVF_NT_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[OFF_HDF_IVF_NT_CM], file = "OFF_HDF_IVF_NT_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[OFF_HDF_IVF_NT_EC], file = "OFF_HDF_IVF_NT_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[OFF_HDF_IVF_iPSC_ESC], file = "OFF_HDF_IVF_iPSC_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[OFF_HDF_IVF_iPSC_CM], file = "OFF_HDF_IVF_iPSC_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[OFF_HDF_IVF_iPSC_EC], file = "OFF_HDF_IVF_iPSC_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[RD_HDF_IVF_NT_ESC], file = "RD_HDF_IVF_NT_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RD_HDF_IVF_NT_CM], file = "RD_HDF_IVF_NT_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RD_HDF_IVF_NT_EC], file = "RD_HDF_IVF_NT_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[RD_HDF_IVF_iPSC_ESC], file = "RD_HDF_IVF_iPSC_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RD_HDF_IVF_iPSC_CM], file = "RD_HDF_IVF_iPSC_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RD_HDF_IVF_iPSC_EC], file = "RD_HDF_IVF_iPSC_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[RU_HDF_IVF_NT_ESC], file = "RU_HDF_IVF_NT_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RU_HDF_IVF_NT_CM], file = "RU_HDF_IVF_NT_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RU_HDF_IVF_NT_EC], file = "RU_HDF_IVF_NT_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(rownames(res21)[RU_HDF_IVF_iPSC_ESC], file = "RU_HDF_IVF_iPSC_ESC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RU_HDF_IVF_iPSC_CM], file = "RU_HDF_IVF_iPSC_CM.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res21)[RU_HDF_IVF_iPSC_EC], file = "RU_HDF_IVF_iPSC_EC.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)


write.table(fpms, file = "FPM_ext.csv",quote=FALSE,sep="\t")
write.table(fpkms, file = "FPKM_ext.csv",quote=FALSE,sep="\t")





D1 <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res9$log2FoldChange,res9$padj))
p <- ggplot(D1, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT_CM_vs_IVF.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res13$log2FoldChange,res13$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("EC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT_EC_vs_IVF.pdf")
dev.off()

D <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res17$log2FoldChange,res17$padj))
p <- ggplot(D, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("ESC log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
ggsave("NT_ESC_vs_IVF.pdf")
dev.off()




library("GenomicFeatures")
refgenes <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
transcripts <- genes(refgenes, columns=c("tx_id", "tx_name"))
#transcripts <- genes(refgenes, columns=c("gene_name"))
tss <- resize(transcripts, width=1, fix='start')


write.table(as.data.frame(tss), file = "TSS.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(as.data.frame(transcripts), file = "gene_body.bed",quote=FALSE,sep="\t",col.names=FALSE)


#resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata1)[1] <- 'gene'
#write.csv(resdata1, file = "HDF_IVF1_CM.csv")

#resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata1)[2] <- 'gene'
#write.csv(resdata2, file = "HDF_IVF1_EC.csv")

#resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata3)[2] <- 'gene'
#write.csv(resdata3, file = "HDF_IVF2_CM.csv")

#resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata4)[2] <- 'gene'
#write.csv(resdata4, file = "HDF_IVF2_EC.csv")

#resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata5)[2] <- 'gene'
#write.csv(resdata5, file = "HDF_NT1_CM.csv")

#resdata6 <- merge(as.data.frame(res6), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata6)[2] <- 'gene'
#write.csv(resdata6, file = "HDF_NT1_EC.csv")

#resdata7 <- merge(as.data.frame(res7), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata7)[2] <- 'gene'
#write.csv(resdata7, file = "HDF_NT2_CM.csv")

#resdata8 <- merge(as.data.frame(res8), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata8)[2] <- 'gene'
#write.csv(resdata8, file = "HDF_IVF1_EC.csv")













#resdata9 <- merge(as.data.frame(res9), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata9)[2] <- 'gene'
#write.csv(resdata9, file = "IVF1_CM_NT1_CM.csv")

#resdata10 <- merge(as.data.frame(res10), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata10)[2] <- 'gene'
#write.csv(resdata10, file = "IVF2_CM_NT1_CM.csv")

#resdata11 <- merge(as.data.frame(res11), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata11)[2] <- 'gene'
#write.csv(resdata10, file = "IVF1_CM_NT2_CM.csv")

#resdata12 <- merge(as.data.frame(res12), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata12)[2] <- 'gene'
#write.csv(resdata10, file = "IVF2_CM_NT2_CM.csv")


#resdata13 <- merge(as.data.frame(res13), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata13)[2] <- 'gene'
#write.csv(resdata13, file = "IVF1_CM_NT1_EC.csv")

#resdata114 <- merge(as.data.frame(res14), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata14)[2] <- 'gene'
#write.csv(resdata14, file = "IVF2_CM_NT1_EC.csv")

#resdata15 <- merge(as.data.frame(res15), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata15)[2] <- 'gene'
#write.csv(resdata15, file = "IVF1_CM_NT2_EC.csv")

#resdata16 <- merge(as.data.frame(res16), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
#names(resdata16)[2] <- 'gene'
#write.csv(resdata16, file = "IVF2_CM_NT2_EC.csv")



#Normalised counts to tab delimited file for GSEA, etc.
#write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
#mcols(res, use.names = T)
#write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

#Replacing outlier value with estimated value as predicted by distrubution using "trimmed mean" approach. 
#ddsClean <- replaceOutliersWithTrimmedMean(dds)
#ddsClean <- DESeq(ddsClean)
#tab <- table(initial = results(dds)$padj < 0.05,cleaned = results(ddsClean)$padj < 0.05)
#addmargins(tab)
#write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
#resClean <- results(ddsClean)
#resClean = subset(res, padj<0.05)
#resClean <- resClean[order(resClean$padj),]
#write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
#plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
#dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
#dev.off()

#Transform raw counts into normalized values. DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
#Variance stabilization is very good for heatmaps, etc.
#rld <- rlogTransformation(dds, blind=F)
#vsd <- varianceStabilizingTransformation(dds, blind=F)

# save normalized values
#write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts-replicate1.txt"), sep = '\t'))
#write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts-replicate1.txt"), sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
#par(mai = ifelse(1:10 <= 2, par('mai'),0))
#px <- counts(dds)[,1] / sizeFactors(dds)[1]
#ord <- order(px)
#ord <- ord[px[ord] < 150]
#ord <- ord[seq(1,length(ord),length=50)]
#last <- ord[length(ord)]
#vstcol <- c('blue','black')
#matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
#legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
#dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
#dev.off()

#Do some clustering ... use transfored datasets
#library("RColorBrewer")
#library("gplots")
#distsRL <- dist(t(assay(rld)))
#mat <- as.matrix(distsRL)
#rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
#dev.copy(png, paste0(outputPrefix, "-rld-clustering1-replicate1.png"))
#dev.off()

#Do some clustering ...
#distsRL <- dist(t(assay(vsd2)))
#mat <- as.matrix(distsRL)
#rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
##Or if you want conditions use:
##rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
#dev.copy(png, paste0(outputPrefix, "-vsd-clustering2-replicate1.png"))
#dev.off()

#Now to PCA ...
#library("genefilter")
#library("ggplot2")
#library("grDevices")


#rv <- rowVars(assay(vsd))
#select <- order(rv, decreasing=T)[seq_len(min(10000,length(rv)))]
#pc <- prcomp(t(assay(vsd)[select,]))
## set condition
#condition <- sampleCondition
#scores <- data.frame(pc$x, condition)
#(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
#+ geom_point(size = 5)
#+ ggtitle("Principal Components")
#+ scale_colour_brewer(name = " ", palette = "Set1")
#+ theme(
#  plot.title = element_text(face = 'bold'),
#  legend.position = c(.9,.2),
#  legend.key = element_rect(fill = 'NA'),
#  legend.text = element_text(size = 10, face = "bold"),
#  axis.text.y = element_text(colour = "Black"),
#  axis.text.x = element_text(colour = "Black"),
#  axis.title.x = element_text(face = "bold"),
#  axis.title.y = element_text(face = 'bold'),
#  panel.grid.major.x = element_blank(),
#  panel.grid.major.y = element_blank(),
#  panel.grid.minor.x = element_blank(),
#  panel.grid.minor.y = element_blank(),
#  panel.background = element_rect(color = 'black',fill = NA)
#))
#ggsave(pcaplot,file=paste0(outputPrefix, "-vsd-ggplot2-replicate1.pdf"))
#
#pcaVar1 <- round(100*pc$sdev^2/sum(pc$sdev^2),1)


#require(gridExtra)
# set condition
#condition <- sampleCondition
#scores <- data.frame(pc$x, condition)
#(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
#+ geom_point(size = 5)
#+ ggtitle("Principal Components")
#+ scale_colour_brewer(name = " ", palette = "Set1")
#+ theme(
#  plot.title = element_text(face = 'bold'),
#  legend.position = c(.9,.2),
#  legend.key = element_rect(fill = 'NA'),
#  legend.text = element_text(size = 10, face = "bold"),
#  axis.text.y = element_text(colour = "Black"),
#  axis.text.x = element_text(colour = "Black"),
#  axis.title.x = element_text(face = "bold"),
#  axis.title.y = element_text(face = 'bold'),
#  panel.grid.major.x = element_blank(),
#  panel.grid.major.y = element_blank(),
#  panel.grid.minor.x = element_blank(),
#  panel.grid.minor.y = element_blank(),
#  panel.background = element_rect(color = 'black',fill = NA)
#))

## set condition
#condition <- sampleCondition
#scores <- data.frame(pc$x, condition)
#(pcaplot2 <- ggplot(scores, aes(x = PC1, y = PC3, col = (factor(condition))))
#+ geom_point(size = 5)
#+ ggtitle("Principal Components")
#+ scale_colour_brewer(name = " ", palette = "Set1")
#+ theme(
#  plot.title = element_text(face = 'bold'),
#  legend.position = c(.9,.2),
#  legend.key = element_rect(fill = 'NA'),
#  legend.text = element_text(size = 10, face = "bold"),
#  axis.text.y = element_text(colour = "Black"),
#  axis.text.x = element_text(colour = "Black"),
#  axis.title.x = element_text(face = "bold"),
#  axis.title.y = element_text(face = 'bold'),
#  panel.grid.major.x = element_blank(),
#  panel.grid.major.y = element_blank(),
#  panel.grid.minor.x = element_blank(),
#  panel.grid.minor.y = element_blank(),
#  panel.background = element_rect(color = 'black',fill = NA)
#))



#rv <- rowVars(assay(rld))
#select <- order(rv, decreasing=T)[seq_len(min(2000,length(rv)))]
#pc <- prcomp(t(assay(rld)[select,]))
## set condition
#condition <- sampleCondition
#scores <- data.frame(pc$x, condition)
#(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
#+ geom_point(size = 5)
#+ ggtitle("Principal Components")
#+ scale_colour_brewer(name = " ", palette = "Set1")
#+ theme(
#  plot.title = element_text(face = 'bold'),
#  legend.position = c(.9,.2),
#  legend.key = element_rect(fill = 'NA'),
#  legend.text = element_text(size = 10, face = "bold"),
#  axis.text.y = element_text(colour = "Black"),
#  axis.text.x = element_text(colour = "Black"),
#  axis.title.x = element_text(face = "bold"),
#  axis.title.y = element_text(face = 'bold'),
#  panel.grid.major.x = element_blank(),
#  panel.grid.major.y = element_blank(),
#  panel.grid.minor.x = element_blank(),
#  panel.grid.minor.y = element_blank(),
#  panel.background = element_rect(color = 'black',fill = NA)
#))
#ggsave(pcaplot,file=paste0(outputPrefix, "-rld-ggplot2-replicate1.pdf"))

#pcaVar2 <- round(100*pc$sdev^2/sum(pc$sdev^2),1)

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
#head(assay(rld))
#par(mfrow=c(1,3))
#plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='ES')
#plot(log2(1+counts(dds,normalized=T)[,3:4]),col='black',pch=20,cex=0.3, main='DMSO')
#plot(log2(1+counts(dds,normalized=T)[,5:6]),col='black',pch=20,cex=0.3, main='iBET')
##plot(log2(1+counts(dds,normalized=T)[,7:8]),col='black',pch=20,cex=0.3, main='d3 DMSO')
##plot(log2(1+counts(dds,normalized=T)[,9:10]),col='black',pch=20,cex=0.3, main='d3 iBET')
#dev.copy(png,"counts-reps-comp.png")
#dev.off()

#par(mfrow=c(1,3))
#plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "ES")
#plot(assay(rld)[,3:4],col='#00000020',pch=20,cex=0.3, main = "DMSO")
#plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "iBET")
##plot(assay(rld)[,7:8],col='#00000020',pch=20,cex=0.3, main = "d3 DMSO")
##plot(assay(rld)[,9:10],col='#00000020',pch=20,cex=0.3, main = "d3 iBET")
#dev.copy(png,"norm-reps-comp.png") 
#dev.off()

