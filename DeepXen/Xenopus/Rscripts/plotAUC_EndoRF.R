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
s
fpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/fpr_kears6_TL.npy")
tpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/tpr_kears6_TL.npy")
s
#Get as a data frame
F1 <- data.frame(FPR = fpr_full, TPR = tpr_full)
F1$ID <- "Full_seed1"

F2 <- data.frame(FPR = fpr_neg, TPR = tpr_neg)
F2$ID <- "Neg_seed1"

F3 <- data.frame(FPR = fpr_pos, TPR = tpr_pos)
F3$ID <- "Pos_seed1"

F4 <- data.frame(FPR = fpr_full_2, TPR = tpr_full_2)
F4$ID <- "Full_seed2"

F5 <- data.frame(FPR = fpr_neg_2, TPR = tpr_neg_2)
F5$ID <- "Neg_seed2"

F6 <- data.frame(FPR = fpr_pos_2, TPR = tpr_pos_2)
F6$ID <- "Pos_seed2"

F7 <- data.frame(FPR = fpr_tlb, TPR = tpr_tlb)
F7$ID <- "TLB_seed1"

F8 <- data.frame(FPR = fpr_tlb_2, TPR = tpr_tlb_2)
F8$ID <- "TLB_seed2"

F9 <- data.frame(FPR = fpr_tl, TPR = tpr_tl)
F9$ID <- "TL_seed1"


F10 <- data.frame(FPR = fpr_tl_2, TPR = tpr_tl_2)
F10$ID <- "TL_seed2"

#Merge data frame
AUC <- rbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10)

#Do plotting
library(ggplot2)
p <- ggplot(AUC, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("dashed","dashed","solid","solid","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

ggsave(filename=paste("/Volumes/Overflow/EndoFull/AUC.pdf",sep=""),width = 20, height = 20, plot = p)






fpr_full <- np$load("/Volumes/Overflow/EndoFull/pr_kearst.npy")
tpr_full <- np$load("/Volumes/Overflow/EndoFull/re_kearst.npy")

fpr_neg <- np$load("/Volumes/Overflow/EndoFull/pr_kears1.npy")
tpr_neg <- np$load("/Volumes/Overflow/EndoFull/re_kears1.npy")

fpr_pos <- np$load("/Volumes/Overflow/EndoFull/pr_kears2.npy")
tpr_pos <- np$load("/Volumes/Overflow/EndoFull/re_kears2.npy")

fpr_full_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears4.npy")
tpr_full_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears4.npy")

fpr_neg_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears5.npy")
tpr_neg_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears5.npy")

fpr_pos_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears6.npy")
tpr_pos_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears6.npy")

fpr_tlb <- np$load("/Volumes/Overflow/EndoFull/pr_kears_TL.npy")
tpr_tlb <- np$load("/Volumes/Overflow/EndoFull/re_kears_TL.npy")

fpr_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears1_TL.npy")
tpr_tlb_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears1_TL.npy")

fpr_tl <- np$load("/Volumes/Overflow/EndoFull/pr_kears5_TL.npy")
tpr_tl <- np$load("/Volumes/Overflow/EndoFull/re_kears5_TL.npy")

fpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/pr_kears6_TL.npy")
tpr_tl_2 <- np$load("/Volumes/Overflow/EndoFull/re_kears6_TL.npy")

#Get as a data frame
F1 <- data.frame(FPR = fpr_full, TPR = tpr_full)
F1$ID <- "Full_seed1"

F2 <- data.frame(FPR = fpr_neg, TPR = tpr_neg)
F2$ID <- "Neg_seed1"

F3 <- data.frame(FPR = fpr_pos, TPR = tpr_pos)
F3$ID <- "Pos_seed1"

F4 <- data.frame(FPR = fpr_full_2, TPR = tpr_full_2)
F4$ID <- "Full_seed2"

F5 <- data.frame(FPR = fpr_neg_2, TPR = tpr_neg_2)
F5$ID <- "Neg_seed2"

F6 <- data.frame(FPR = fpr_pos_2, TPR = tpr_pos_2)
F6$ID <- "Pos_seed2"

F7 <- data.frame(FPR = fpr_tlb, TPR = tpr_tlb)
F7$ID <- "TLB_seed1"

F8 <- data.frame(FPR = fpr_tlb_2, TPR = tpr_tlb_2)
F8$ID <- "TLB_seed2"

F9 <- data.frame(FPR = fpr_tl, TPR = tpr_tl)
F9$ID <- "TL_seed1"


F10 <- data.frame(FPR = fpr_tl_2, TPR = tpr_tl_2)
F10$ID <- "TL_seed2"

#Merge data frame
AUC <- rbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10)

#Do plotting
library(ggplot2)
p <- ggplot(AUC, aes(x=TPR, y=FPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("dashed","dashed","solid","solid","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

ggsave(filename=paste("/Volumes/Overflow/EndoFull/PR.pdf",sep=""),width = 20, height = 20, plot = p)

