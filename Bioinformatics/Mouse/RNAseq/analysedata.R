
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")
library(dplyr)
library(ggplot2)
library("ggbeeswarm")
library("RColorBrewer")
library("gplots")
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
geneLength <- as.data.frame(unlist(exonic.gene.sizes))



sampleNames <- c("MEF1","MEF2","MEF3","1C1nt","1C2nt","1C3nt","2C1nt","2C2nt","2C3nt","4C1nt","4C2nt","4C3nt","8C1nt","8C2nt","8C3nt","M1nt","M2nt","M3nt","B1nt","B2nt","B3nt","1C1","1C2","1C3","2C1","2C2","2C3","4C1","4C2","4C3","8C1","8C2","8C3","M1","M2","M3","B1","B2","B3")
sampleFiles <- c("SRR7007704/SRR7007704.counts","SRR7007705/SRR7007705.counts","SRR7007706/SRR7007706.counts","SRR7007725/SRR7007725.counts","SRR7007726/SRR7007726.counts","SRR7007727/SRR7007727.counts","SRR7007728/SRR7007728.counts","SRR7007729/SRR7007729.counts","SRR7007730/SRR7007730.counts","SRR7007731/SRR7007731.counts","SRR7007732/SRR7007732.counts","SRR7007733/SRR7007733.counts","SRR7007734/SRR7007734.counts","SRR7007735/SRR7007735.counts","SRR7007736/SRR7007736.counts","SRR7007740/SRR7007740.counts","SRR7007741/SRR7007741.counts","SRR7007742/SRR7007742.counts","SRR7007737/SRR7007737.counts","SRR7007738/SRR7007738.counts","SRR7007739/SRR7007739.counts","SRR7007686/SRR7007686.counts","SRR7007687/SRR7007687.counts","SRR7007688/SRR7007688.counts","SRR7007689/SRR7007689.counts","SRR7007690/SRR7007690.counts","SRR7007691/SRR7007691.counts","SRR7007692/SRR7007692.counts","SRR7007693/SRR7007693.counts","SRR7007694/SRR7007694.counts","SRR7007695/SRR7007695.counts","SRR7007696/SRR7007696.counts","SRR7007697/SRR7007697.counts","SRR7007701/SRR7007701.counts","SRR7007702/SRR7007702.counts","SRR7007703/SRR7007703.counts","SRR7007698/SRR7007698.counts","SRR7007699/SRR7007699.counts","SRR7007700/SRR7007700.counts")
sampleCondition <- c("MEF","MEF","MEF","NT1C","NT1C","NT1C","NT2C","NT2C","NT2C","NT4C","NT4C","NT4C","NT8C","NT8C","NT8C","NTM","NTM","NTM","NTB","NTB","NTB","ivf1C","ivf1C","ivf1C","ivf2C","ivf2C","ivf2C","ivf4C","ivf4C","ivf4C","ivf8C","ivf8C","ivf8C","ivfM","ivfM","ivfM","ivfB","ivfB","ivfB")
sampleReplicate <- c("1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3")
#sampleNames <- sampleNames[1:8]
#sampleFiles <- sampleFiles[1:8]
#sampleCondition <-  sampleCondition[1:8]
#sampleReplicate <-sampleReplicate[1:8]
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, rep = sampleReplicate)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/mnt/scratch/gurdon/cap76/DeepXen/Mouse/RNAseq/", design=~condition)


#1 cell
#SRR7007725 SRR7007726 SRR7007727
#2
#SRR7007728 SRR7007729 SRR7007730
#4
#SRR7007731 SRR7007732 SRR7007733
#8
#SRR7007734 SRR7007735 SRR7007736
#M
#SRR7007740 SRR7007741 SRR7007742
#B
#SRR7007737 SRR7007738 SRR7007739

#SRR7007686 SRR7007687 SRR7007688
#SRR7007689 SRR7007690 SRR7007691
#SRR7007692 SRR7007693 SRR7007694
#SRR7007695 SRR7007696 SRR7007697
#SRR7007701 SRR7007702 SRR7007703
#SRR7007698 SRR7007699 SRR7007700

directory <- "."
setwd(directory)


#sampleNames <- c("MEF1","MEF2","MEF3","1C1nt1","1C1nt2","1C1nt3") #,"2C1nt","2C2nt","2C3nt","4C1nt","4C2nt","4C3nt","8C1nt","8C2nt","8C3nt","M1nt","M2nt","M3nt","B1nt","B2nt","B3nt","1C1","1C1","1C1","2C1","2C2","2C3","4C1","4C2","4C3","8C1","8C2","8C3","M1","M2","M3","B1","B2","B3")
#sampleFiles <- c("SRR7007686/SRR7007686.counts","SRR7007705/SRR7007705.counts","SRR7007706/SRR7007706.counts","SRR7007725/SRR7007725.counts","SRR7007726/SRR7007726.counts","SRR7007727/SRR7007727.counts") #,"SRR7007728/SRR7007728.counts","SRR7007729/SRR7007729.counts","SRR7007730/SRR7007730.counts","SRR7007731/SRR7007731.counts","SRR7007732/SRR7007732.counts","SRR7007733/SRR7007733.counts","SRR7007734/SRR7007734.counts","SRR7007735/SRR7007735.counts","SRR7007736/SRR7007736.counts","SRR7007740/SRR7007740.counts","SRR7007741/SRR7007741.counts","SRR7007742/SRR7007742.counts","SRR7007737/SRR7007737.counts","SRR7007738/SRR7007738.counts","SRR7007739/SRR7007739.counts","SRR7007686/SRR7007686.counts","SRR7007687/SRR7007687.counts","SRR7007688/SRR7007688.counts","SRR7007689/SRR7007689.counts","SRR7007690/SRR7007690.counts","SRR7007691/SRR7007691.counts","SRR7007692/SRR7007692.counts","SRR7007693/SRR7007693.counts","SRR7007694/SRR7007694.counts","SRR7007695/SRR7007695.counts","SRR7007696/SRR7007696.counts","SRR7007697/SRR7007697.counts","SRR7007701/SRR7007701.counts","SRR7007702/SRR7007702.counts","SRR7007703/SRR7007703.counts","SRR7007698/SRR7007698.counts","SRR7007699/SRR7007699.counts","SRR7007700/SRR7007700.counts")
#sampleCondition <- c("MEF","MEF","MEF","NT1C","NT1C","NT1C") #,"NT2C","NT2C","NT2C","NT4C","NT4C","NT4C","NT8C","NT8C","NT8C","NTM","NTM","NTM","NTB","NTB","NTB","ivf1C","ivf1C","ivf1C","ivf2C","ivf2C","ivf2C","ivf4C","ivf4C","ivf4C","ivf8C","ivf8C","ivf8C","ivfM","ivfM","ivfM","ivfB","ivfB","ivfB")
#sampleReplicate <- c("1","2","3","1","2","3") #,"1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3","1","2","3")
#sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, rep = sampleReplicate)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/mnt/scratch/gurdon/cap76/DeepXen/Mouse/RNAseq/", design=~condition)

#Now we get gene lengths
#txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
#exons.list.per.gene <- exonsBy(txdb,by="gene")
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
#geneLength <- as.data.frame(unlist(exonic.gene.sizes))

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
res1<-results(dds, contrast=c("condition","MEF","ivf1C"))
res2<-results(dds, contrast=c("condition","MEF","ivf2C"))
res3<-results(dds, contrast=c("condition","MEF","ivf4C"))
res4<-results(dds, contrast=c("condition","MEF","ivf8C"))
res5<-results(dds, contrast=c("condition","MEF","ivfM"))
res6<-results(dds, contrast=c("condition","MEF","ivfB"))

res7<-results(dds, contrast=c("condition","MEF","NT1C"))
res8<-results(dds, contrast=c("condition","MEF","NT2C"))
res9<-results(dds, contrast=c("condition","MEF","NT4C"))
res10<-results(dds, contrast=c("condition","MEF","NT8C"))
res11<-results(dds, contrast=c("condition","MEF","NTM"))
res12<-results(dds, contrast=c("condition","MEF","NTB"))

res13<-results(dds, contrast=c("condition","NT1C","ivf1C"))
res14<-results(dds, contrast=c("condition","NT2C","ivf2C"))
res15<-results(dds, contrast=c("condition","NT4C","ivf4C"))
res16<-results(dds, contrast=c("condition","NT8C","ivf8C"))
res17<-results(dds, contrast=c("condition","NTM","ivfM"))
res18<-results(dds, contrast=c("condition","NTB","ivfB"))

ON_MEF_IVF_NT_1C <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_1C <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res13)$pval<0.05 & as.data.frame(res13)$log2FoldChange<0)
RD_MEF_IVF_NT_1C <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange>0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_1C <- which(as.data.frame(res1)$pval<0.05 & as.data.frame(res1)$log2FoldChange<0 & as.data.frame(res7)$pval<0.05 & as.data.frame(res7)$log2FoldChange<0 & as.data.frame(res13)$pval>0.05)

ON_MEF_IVF_NT_2C <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_2C <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res14)$pval<0.05 & as.data.frame(res14)$log2FoldChange<0)
RD_MEF_IVF_NT_2C <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange>0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_2C <- which(as.data.frame(res2)$pval<0.05 & as.data.frame(res2)$log2FoldChange<0 & as.data.frame(res8)$pval<0.05 & as.data.frame(res8)$log2FoldChange<0 & as.data.frame(res14)$pval>0.05)

ON_MEF_IVF_NT_4C <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res15)$pval<0.05 & as.data.frame(res15)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_4C <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res15)$pval<0.05 & as.data.frame(res15)$log2FoldChange<0)
RD_MEF_IVF_NT_4C <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange>0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_4C <- which(as.data.frame(res3)$pval<0.05 & as.data.frame(res3)$log2FoldChange<0 & as.data.frame(res9)$pval<0.05 & as.data.frame(res9)$log2FoldChange<0 & as.data.frame(res15)$pval>0.05)

ON_MEF_IVF_NT_8C <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res16)$pval<0.05 & as.data.frame(res16)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_8C <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res16)$pval<0.05 & as.data.frame(res16)$log2FoldChange<0)
RD_MEF_IVF_NT_8C <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange>0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_8C <- which(as.data.frame(res4)$pval<0.05 & as.data.frame(res4)$log2FoldChange<0 & as.data.frame(res10)$pval<0.05 & as.data.frame(res10)$log2FoldChange<0 & as.data.frame(res16)$pval>0.05)

ON_MEF_IVF_NT_M <- which(as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange>0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_M <- which(as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange<0 & as.data.frame(res17)$pval<0.05 & as.data.frame(res17)$log2FoldChange<0)
RD_MEF_IVF_NT_M <- which(as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange>0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_M <- which(as.data.frame(res5)$pval<0.05 & as.data.frame(res5)$log2FoldChange<0 & as.data.frame(res11)$pval<0.05 & as.data.frame(res11)$log2FoldChange<0 & as.data.frame(res17)$pval>0.05)

ON_MEF_IVF_NT_B <- which(as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange>0 & as.data.frame(res18)$pval<0.05 & as.data.frame(res18)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
OFF_MEF_IVF_NT_B <- which(as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange<0 & as.data.frame(res18)$pval<0.05 & as.data.frame(res18)$log2FoldChange<0)
RD_MEF_IVF_NT_B <- which(as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange>0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange>0 & (fpkms[,1]+fpkms[,2]+fpkms[,3])/3 > 1)
RU_MEF_IVF_NT_B <- which(as.data.frame(res6)$pval<0.05 & as.data.frame(res6)$log2FoldChange<0 & as.data.frame(res12)$pval<0.05 & as.data.frame(res12)$log2FoldChange<0 & as.data.frame(res18)$pval>0.05)


write.table(rownames(res1)[ON_MEF_IVF_NT_1C], file = "ON_MEF_IVF_NT_1C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[ON_MEF_IVF_NT_2C], file = "ON_MEF_IVF_NT_2C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[ON_MEF_IVF_NT_4C], file = "ON_MEF_IVF_NT_4C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[ON_MEF_IVF_NT_8C], file = "ON_MEF_IVF_NT_8C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[ON_MEF_IVF_NT_M], file = "ON_MEF_IVF_NT_M.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[ON_MEF_IVF_NT_B], file = "ON_MEF_IVF_NT_B.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)


write.table(rownames(res1)[OFF_MEF_IVF_NT_1C], file = "OFF_MEF_IVF_NT_1C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[OFF_MEF_IVF_NT_2C], file = "OFF_MEF_IVF_NT_2C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[OFF_MEF_IVF_NT_4C], file = "OFF_MEF_IVF_NT_4C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[OFF_MEF_IVF_NT_8C], file = "OFF_MEF_IVF_NT_8C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[OFF_MEF_IVF_NT_M], file = "OFF_MEF_IVF_NT_M.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[OFF_MEF_IVF_NT_B], file = "OFF_MEF_IVF_NT_B.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)


write.table(rownames(res1)[RU_MEF_IVF_NT_1C], file = "RU_MEF_IVF_NT_1C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RU_MEF_IVF_NT_2C], file = "RU_MEF_IVF_NT_2C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RU_MEF_IVF_NT_4C], file = "RU_MEF_IVF_NT_4C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RU_MEF_IVF_NT_8C], file = "RU_MEF_IVF_NT_8C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RU_MEF_IVF_NT_M], file = "RU_MEF_IVF_NT_M.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RU_MEF_IVF_NT_B], file = "RU_MEF_IVF_NT_B.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)


write.table(rownames(res1)[RD_MEF_IVF_NT_1C], file = "RD_MEF_IVF_NT_1C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RD_MEF_IVF_NT_2C], file = "RD_MEF_IVF_NT_2C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RD_MEF_IVF_NT_4C], file = "RD_MEF_IVF_NT_4C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RD_MEF_IVF_NT_8C], file = "RD_MEF_IVF_NT_8C.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RD_MEF_IVF_NT_M], file = "RD_MEF_IVF_NT_M.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(res1)[RD_MEF_IVF_NT_B], file = "RD_MEF_IVF_NT_B.csv",quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(fpms, file = "FPM_ext.csv",quote=FALSE,sep="\t")
write.table(fpkms, file = "FPKM_ext.csv",quote=FALSE,sep="\t")


#D1 <- as.data.frame(cbind( log2( (fpkms[,1]+fpkms[,2]+fpkms[,3])/3  +1 )  ,res9$log2FoldChange,res9$padj))
#p <- ggplot(D1, aes_string(x="V1", y="V2")) + theme_classic() + geom_quasirandom(alpha = 0.1, width = 0.3) + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") + xlab("HDF log2(FPKM +1)") + ylab("CM log2FC(NT vs IVF)") + geom_hline(yintercept=c(-1,1),linetype="dashed", color = "red")
#ggsave("NT_CM_vs_IVF.pdf")
#dev.off()

#write.table(as.data.frame(tss), file = "TSS.bed",quote=FALSE,sep="\t",col.names=FALSE)
#write.table(as.data.frame(transcripts), file = "gene_body.bed",quote=FALSE,sep="\t",col.names=FALSE)


library("GenomicFeatures")
refgenes <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf",format="gtf")
transcripts <- genes(refgenes, columns=c("tx_id", "tx_name"))
#transcripts <- genes(refgenes, columns=c("gene_name"))
tss <- resize(transcripts, width=1, fix='start')

write.table(as.data.frame(tss), file = "TSS.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(as.data.frame(transcripts), file = "gene_body.bed",quote=FALSE,sep="\t",col.names=FALSE)

