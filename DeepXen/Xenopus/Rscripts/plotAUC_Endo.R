library(reticulate)
library(pheatmap)

np <- import("numpy")

fpr_full <- np$load("/Volumes/Overflow/EndoFull/fpr_kearst.npy")
tpr_full <- np$load("/Volumes/Overflow/EndoFull/tpr_kearst.npy")

fpr_neg <- np$load("/Volumes/Overflow/EndoFull/fpr_kears1.npy")
tpr_neg <- np$load("/Volumes/Overflow/EndoFull/tpr_kears1.npy")

fpr_pos <- np$load("/Volumes/Overflow/EndoFull/fpr_kears2.npy")
tpr_pos <- np$load("/Volumes/Overflow/EndoFull/tpr_kears2.npy")

fpr_full_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears4.npy")
tpr_full_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears4.npy")

fpr_neg_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears5.npy")
tpr_neg_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears5.npy")

fpr_pos_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears6.npy")
tpr_pos_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears6.npy")

fpr_tlb <- np$load("/Volumes/Overflow/EndoFull/fpr_kears_TL.npy")
tpr_tlb <- np$load("/Volumes/Overflow/EndoFull/tpr_kears_TL.npy")

fpr_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears1_TL.npy")
tpr_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears1_TL.npy")

fpr_tl <- np$load("/Volumes/Overflow/EndoFull/fpr_kears5_TL.npy")
tpr_tl <- np$load("/Volumes/Overflow/EndoFull/tpr_kears5_TL.npy")

fpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears6_TL.npy")
tpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears6_TL.npy")

#Get as a data frame
F1 <- data.frame(FPR = fpr_full, TPR = tpr_full)
F1$ID <- "Combined (seed 1)"

F2 <- data.frame(FPR = fpr_neg, TPR = tpr_neg)
F2$ID <- "Chromatin (seed 1)"

F3 <- data.frame(FPR = fpr_pos, TPR = tpr_pos)
F3$ID <- "Expression (seed 1)"

F4 <- data.frame(FPR = fpr_full_2, TPR = tpr_full_2)
F4$ID <- "Combined (seed 2)"

F5 <- data.frame(FPR = fpr_neg_2, TPR = tpr_neg_2)
F5$ID <- "Chromatin (seed 2)"

F6 <- data.frame(FPR = fpr_pos_2, TPR = tpr_pos_2)
F6$ID <- "Expression (seed 2)"

F7 <- data.frame(FPR = fpr_tlb, TPR = tpr_tlb)
F7$ID <- "Naive TL (seed 1)"

F8 <- data.frame(FPR = fpr_tlb_2, TPR = tpr_tlb_2)
F8$ID <- "Naive TL (seed 2)"

F9 <- data.frame(FPR = fpr_tl, TPR = tpr_tl)
F9$ID <- "TL (seed 1)"

F10 <- data.frame(FPR = fpr_tl_2, TPR = tpr_tl_2)
F10$ID <- "TL (seed 2)"

#Merge data frame
AUC <- rbind(F1,F2,F3,F4,F5,F6)

#Do plotting
library(ggplot2)
p <- ggplot(AUC[,], aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","dashed","solid","dashed","solid","dashed","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#88419d","#4d004b","#0570b0","#023858","#d7301f","#7f0000")) #+ geom_line(size=2)

ggsave(filename=paste("/Volumes/Overflow/EndoFull/AUC.pdf",sep=""),width = 10, height = 10, plot = p)

pr_full <- np$load("/Volumes/Overflow/EndoFull/pr_kearst.npy")
re_full <- np$load("/Volumes/Overflow/EndoFull/re_kearst.npy")

pr_neg <- np$load("/Volumes/Overflow/EndoFull/pr_kears1.npy")
re_neg <- np$load("/Volumes/Overflow/EndoFull/re_kears1.npy")

pr_pos <- np$load("/Volumes/Overflow/EndoFull/pr_kears2.npy")
re_pos <- np$load("/Volumes/Overflow/EndoFull/re_kears2.npy")

pr_full_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears4.npy")
re_full_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears4.npy")

pr_neg_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears5.npy")
re_neg_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears5.npy")

pr_pos_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears6.npy")
re_pos_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears6.npy")

pr_tlb <- np$load("/Volumes/Overflow/EndoFull/pr_kears_TL.npy")
re_tlb <- np$load("/Volumes/Overflow/EndoFull/re_kears_TL.npy")

pr_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears1_TL.npy")
re_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears1_TL.npy")

pr_tl <- np$load("/Volumes/Overflow/EndoFull/pr_kears5_TL.npy")
re_tl <- np$load("/Volumes/Overflow/EndoFull/re_kears5_TL.npy")

pr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears6_TL.npy")
re_tl_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears6_TL.npy")

#Get as a data frame
F1 <- data.frame(Pr = pr_full, Re = re_full)
F1$ID <- "Combined (seed 1)"

F2 <- data.frame(Pr = pr_neg, Re = re_neg)
F2$ID <- "Chromatin (seed 1)"

F3 <- data.frame(Pr = pr_pos, Re = re_pos)
F3$ID <- "Expression (seed 1)"

F4 <- data.frame(Pr = pr_full_2, Re = re_full_2)
F4$ID <- "Combined (seed 2)"

F5 <- data.frame(Pr = pr_neg_2, Re  = re_neg_2)
F5$ID <- "Chromatin (seed 3)"

F6 <- data.frame(Pr = pr_pos_2, Re = re_pos_2)
F6$ID <- "Expression (seed 2)"

F7 <- data.frame(Pr = pr_tlb, Re = re_tlb)
F7$ID <- "TLB_seed1"

F8 <- data.frame(Pr = pr_tlb_2, Re = re_tlb_2)
F8$ID <- "TLB_seed2"

F9 <- data.frame(Pr = pr_tl, Re = re_tl)
F9$ID <- "TL_seed1"


F10 <- data.frame(Pr = pr_tl_2, Re = re_tl_2)
F10$ID <- "TL_seed2"

#Merge data frame
PrRe <- rbind(F1,F2,F3,F4,F5,F6) ##,F7,F8,F9,F10)

#Do plotting
#library(ggplot2)
p <- ggplot(PrRe, aes(x=Re, y=Pr, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","dashed","solid","dashed","solid","dashed","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#88419d","#4d004b","#0570b0","#023858","#d7301f","#7f0000")) #+ geom_line(size=2)

ggsave(filename=paste("/Volumes/Overflow/EndoFull/PR.pdf",sep=""),width = 10, height = 10, plot = p)
#+ geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("dashed","dashed","solid","solid","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#ggsave(filename=paste("/Volumes/Overflow/EndoFull/PR.pdf",sep=""),width = 20, height = 20, plot = p)

