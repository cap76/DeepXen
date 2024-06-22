library(reticulate)
library(pheatmap)
library(ggplot2)
require(gridExtra)
library("pracma")

rm(list = ls())
np <- import("numpy")

#First load the ROC curves for the random forest extended control (chrom only)

#First load the ROC curves for the random forest (full model)
rf_fpr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears0_1.npy")
rf_fpr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears0_2.npy")
rf_fpr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears0_3.npy")
rf_fpr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears0_4.npy")
rf_fpr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears0_5.npy")
rf_fpr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears1_1.npy")
rf_fpr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears1_2.npy")
rf_fpr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears1_3.npy")
rf_fpr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears1_4.npy")
rf_fpr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_fpr_kears1_5.npy")

rf_tpr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears0_1.npy")
rf_tpr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears0_2.npy")
rf_tpr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears0_3.npy")
rf_tpr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears0_4.npy")
rf_tpr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears0_5.npy")
rf_tpr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears1_1.npy")
rf_tpr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears1_2.npy")
rf_tpr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears1_3.npy")
rf_tpr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears1_4.npy")
rf_tpr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/rf_tpr_kears1_5.npy")

#Now the NNs
fpr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears0_1.npy")
fpr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears0_2.npy")
fpr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears0_3.npy")
fpr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears0_4.npy")
fpr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears0_5.npy")
fpr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears1_1.npy")
fpr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears1_2.npy")
fpr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears1_3.npy")
fpr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears1_4.npy")
fpr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears1_5.npy")
fpr_full_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears2_1.npy")
fpr_full_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears2_2.npy")
fpr_full_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears2_3.npy")
fpr_full_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears2_4.npy")
fpr_full_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears2_5.npy")
fpr_full_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears4_1.npy")
fpr_full_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears4_2.npy")
fpr_full_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears4_3.npy")
fpr_full_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears4_4.npy")
fpr_full_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears4_5.npy")
fpr_full_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears5_1.npy")
fpr_full_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears5_2.npy")
fpr_full_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears5_3.npy")
fpr_full_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears5_4.npy")
fpr_full_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears5_5.npy")
fpr_full_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears6_1.npy")
fpr_full_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears6_2.npy")
fpr_full_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears6_3.npy")
fpr_full_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears6_4.npy")
fpr_full_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/fpr_kears6_5.npy")

tpr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears0_1.npy")
tpr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears0_2.npy")
tpr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears0_3.npy")
tpr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears0_4.npy")
tpr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears0_5.npy")
tpr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears1_1.npy")
tpr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears1_2.npy")
tpr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears1_3.npy")
tpr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears1_4.npy")
tpr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears1_5.npy")
tpr_full_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears2_1.npy")
tpr_full_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears2_2.npy")
tpr_full_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears2_3.npy")
tpr_full_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears2_4.npy")
tpr_full_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears2_5.npy")
tpr_full_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears4_1.npy")
tpr_full_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears4_2.npy")
tpr_full_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears4_3.npy")
tpr_full_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears4_4.npy")
tpr_full_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears4_5.npy")
tpr_full_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears5_1.npy")
tpr_full_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears5_2.npy")
tpr_full_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears5_3.npy")
tpr_full_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears5_4.npy")
tpr_full_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears5_5.npy")
tpr_full_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears6_1.npy")
tpr_full_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears6_2.npy")
tpr_full_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears6_3.npy")
tpr_full_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears6_4.npy")
tpr_full_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/tpr_kears6_5.npy")

#Now the RF for the GB 
rf_fpr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears0_1.npy")
rf_fpr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears0_2.npy")
rf_fpr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears0_3.npy")
rf_fpr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears0_4.npy")
rf_fpr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears0_5.npy")
rf_fpr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears1_1.npy")
rf_fpr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears1_2.npy")
rf_fpr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears1_3.npy")
rf_fpr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears1_4.npy")
rf_fpr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_fpr_kears1_5.npy")

rf_tpr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears0_1.npy")
rf_tpr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears0_2.npy")
rf_tpr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears0_3.npy")
rf_tpr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears0_4.npy")
rf_tpr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears0_5.npy")
rf_tpr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears1_1.npy")
rf_tpr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears1_2.npy")
rf_tpr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears1_3.npy")
rf_tpr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears1_4.npy")
rf_tpr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_tpr_kears1_5.npy")

#And NN for GB
fpr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears0_1.npy")
fpr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears0_2.npy")
fpr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears0_3.npy")
fpr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears0_4.npy")
fpr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears0_5.npy")
fpr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears1_1.npy")
fpr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears1_2.npy")
fpr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears1_3.npy")
fpr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears1_4.npy")
fpr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears1_5.npy")
fpr_gb_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears2_1.npy")
fpr_gb_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears2_2.npy")
fpr_gb_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears2_3.npy")
fpr_gb_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears2_4.npy")
fpr_gb_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears2_5.npy")
fpr_gb_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears4_1.npy")
fpr_gb_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears4_2.npy")
fpr_gb_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears4_3.npy")
fpr_gb_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears4_4.npy")
fpr_gb_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears4_5.npy")
fpr_gb_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears5_1.npy")
fpr_gb_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears5_2.npy")
fpr_gb_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears5_3.npy")
fpr_gb_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears5_4.npy")
fpr_gb_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears5_5.npy")
fpr_gb_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears6_1.npy")
fpr_gb_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears6_2.npy")
fpr_gb_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears6_3.npy")
fpr_gb_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears6_4.npy")
fpr_gb_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/fpr_kears6_5.npy")

tpr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears0_1.npy")
tpr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears0_2.npy")
tpr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears0_3.npy")
tpr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears0_4.npy")
tpr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears0_5.npy")
tpr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears1_1.npy")
tpr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears1_2.npy")
tpr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears1_3.npy")
tpr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears1_4.npy")
tpr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears1_5.npy")
tpr_gb_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears2_1.npy")
tpr_gb_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears2_2.npy")
tpr_gb_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears2_3.npy")
tpr_gb_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears2_4.npy")
tpr_gb_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears2_5.npy")
tpr_gb_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears4_1.npy")
tpr_gb_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears4_2.npy")
tpr_gb_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears4_3.npy")
tpr_gb_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears4_4.npy")
tpr_gb_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears4_5.npy")
tpr_gb_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears5_1.npy")
tpr_gb_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears5_2.npy")
tpr_gb_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears5_3.npy")
tpr_gb_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears5_4.npy")
tpr_gb_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears5_5.npy")
tpr_gb_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears6_1.npy")
tpr_gb_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears6_2.npy")
tpr_gb_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears6_3.npy")
tpr_gb_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears6_4.npy")
tpr_gb_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/tpr_kears6_5.npy")


#Calculate the AUC
#RFs
D1_1 <- data.frame(FPR = rf_fpr_full_0_1, TPR = rf_tpr_full_0_1)
D1_1$ID <- "RF_1"
D2_1 <- data.frame(FPR = rf_fpr_full_1_1, TPR = rf_tpr_full_1_1)
D2_1$ID <- "RF_2"
D1_1a <- data.frame(FPR = rf_fpr_gb_0_1, TPR = rf_tpr_gb_0_1)
D1_1a$ID <- "RF_GB_1"
D2_1a <- data.frame(FPR = rf_fpr_gb_1_1, TPR = rf_tpr_gb_1_1)
D2_1a$ID <- "RF_GB_2"
#Full NNs
D5_1 <- data.frame(FPR = fpr_full_0_1, TPR = tpr_full_0_1)
D5_1$ID <- "Full_1"
D8_1 <- data.frame(FPR = fpr_full_4_1, TPR = tpr_full_4_1)
D8_1$ID <- "Full_2"
D5_1a <- data.frame(FPR = fpr_gb_0_1, TPR = tpr_gb_0_1)
D5_1a$ID <- "Full_GB_1"
D8_1a <- data.frame(FPR = fpr_gb_4_1, TPR = tpr_gb_4_1)
D8_1a$ID <- "Full_GB_2"

AUC1 <- c(
          trapz(D5_1$FPR,D5_1$TPR),
          trapz(D8_1$FPR,D8_1$TPR),
          trapz(D1_1$FPR,D1_1$TPR),
          trapz(D2_1$FPR,D2_1$TPR),
          trapz(D5_1a$FPR,D5_1a$TPR),
          trapz(D8_1a$FPR,D8_1a$TPR),
          trapz(D1_1a$FPR,D1_1a$TPR),
          trapz(D2_1a$FPR,D2_1a$TPR))

AUC1 <- as.data.frame(AUC1)
rownames(AUC1) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")


D1_2 <- data.frame(FPR = rf_fpr_gb_0_2, TPR = rf_tpr_gb_0_2)
D1_2$ID <- "RF_1"
D2_2 <- data.frame(FPR = rf_fpr_gb_1_2, TPR = rf_tpr_gb_1_2)
D2_2$ID <- "RF_2"
D1_2a <- data.frame(FPR = rf_fpr_gb_0_2, TPR = rf_tpr_gb_0_2)
D1_2a$ID <- "RF_GB_1"
D2_2a <- data.frame(FPR = rf_fpr_gb_1_2, TPR = rf_tpr_gb_1_2)
D2_2a$ID <- "RF_GB_2"

D5_2 <- data.frame(FPR = fpr_full_0_2, TPR = tpr_full_0_2)
D5_2$ID <- "Full_1"
D8_2 <- data.frame(FPR = fpr_full_4_2, TPR = tpr_full_4_2)
D8_2$ID <- "Full_2"

D5_2a <- data.frame(FPR = fpr_gb_0_2, TPR = tpr_gb_0_2)
D5_2a$ID <- "Full_GB_1"
D8_2a <- data.frame(FPR = fpr_gb_4_2, TPR = tpr_gb_4_2)
D8_2a$ID <- "Full_GB_2"

AUC2 <- c(
  trapz(D5_2$FPR,D5_2$TPR),
  trapz(D8_2$FPR,D8_2$TPR),
  trapz(D1_2$FPR,D1_2$TPR),
  trapz(D2_2$FPR,D2_2$TPR),
  trapz(D5_2a$FPR,D5_2a$TPR),
  trapz(D8_2a$FPR,D8_2a$TPR),
  trapz(D1_2a$FPR,D1_2a$TPR),
  trapz(D2_2a$FPR,D2_2a$TPR))

AUC2 <- as.data.frame(AUC2)
rownames(AUC2) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")


D1_3 <- data.frame(FPR = rf_fpr_full_0_3, TPR = rf_tpr_full_0_3)
D1_3$ID <- "RF_1"
D2_3 <- data.frame(FPR = rf_fpr_full_1_3, TPR = rf_tpr_full_1_3)
D2_3$ID <- "RF_2"
D1_3a <- data.frame(FPR = rf_fpr_gb_0_3, TPR = rf_tpr_gb_0_3)
D1_3a$ID <- "RF_GB_1"
D2_3a <- data.frame(FPR = rf_fpr_gb_1_3, TPR = rf_tpr_gb_1_3)
D2_3a$ID <- "RF_GB_2"


D5_3 <- data.frame(FPR = fpr_full_0_3, TPR = tpr_full_0_3)
D5_3$ID <- "Full_1"
D8_3 <- data.frame(FPR = fpr_full_4_3, TPR = tpr_full_4_3)
D8_3$ID <- "Full_2"
D5_3a <- data.frame(FPR = fpr_gb_0_3, TPR = tpr_gb_0_3)
D5_3a$ID <- "Full_GB_1"
D8_3a <- data.frame(FPR = fpr_gb_4_3, TPR = tpr_gb_4_3)
D8_3a$ID <- "Full_GB_2"

#Merge data frame
AUC3 <- c(
  trapz(D5_3$FPR,D5_3$TPR),
  trapz(D8_3$FPR,D8_3$TPR),
  trapz(D1_3$FPR,D1_3$TPR),
  trapz(D2_3$FPR,D2_3$TPR),
  trapz(D5_3a$FPR,D5_3a$TPR),
  trapz(D8_3a$FPR,D8_3a$TPR),
  trapz(D1_3a$FPR,D1_3a$TPR),
  trapz(D2_3a$FPR,D2_3a$TPR))

AUC3 <- as.data.frame(AUC3)
rownames(AUC3) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")

D1_4 <- data.frame(FPR = rf_fpr_full_0_4, TPR = rf_tpr_full_0_4)
D1_4$ID <- "RF_1"
D2_4 <- data.frame(FPR = rf_fpr_full_1_4, TPR = rf_tpr_full_1_4)
D2_4$ID <- "RF_2"
D1_4a <- data.frame(FPR = rf_fpr_gb_0_4, TPR = rf_tpr_gb_0_4)
D1_4a$ID <- "RF_GB_1"
D2_4a <- data.frame(FPR = rf_fpr_gb_1_4, TPR = rf_tpr_gb_1_4)
D2_4a$ID <- "RF_GB_2"

D5_4 <- data.frame(FPR = fpr_full_0_4, TPR = tpr_full_0_4)
D5_4$ID <- "Full_1"
D8_4 <- data.frame(FPR = fpr_full_4_4, TPR = tpr_full_4_4)
D8_4$ID <- "Full_2"
D5_4a <- data.frame(FPR = fpr_gb_0_4, TPR = tpr_gb_0_4)
D5_4a$ID <- "Full_GB_1"
D8_4a <- data.frame(FPR = fpr_gb_4_4, TPR = tpr_gb_4_4)
D8_4a$ID <- "Full_GB_2"
#Merge data frame

AUC4 <- c(
  trapz(D5_4$FPR,D5_4$TPR),
  trapz(D8_4$FPR,D8_4$TPR),
  trapz(D1_4$FPR,D1_4$TPR),
  trapz(D2_4$FPR,D2_4$TPR),
  trapz(D5_4a$FPR,D5_4a$TPR),
  trapz(D8_4a$FPR,D8_4a$TPR),
  trapz(D1_4a$FPR,D1_4a$TPR),
  trapz(D2_4a$FPR,D2_4a$TPR))

AUC4 <- as.data.frame(AUC4)
rownames(AUC4) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")


D1_5 <- data.frame(FPR = rf_fpr_full_0_5, TPR = rf_tpr_full_0_5)
D1_5$ID <- "RF_1"
D2_5 <- data.frame(FPR = rf_fpr_full_1_5, TPR = rf_tpr_full_1_5)
D2_5$ID <- "RF_2"
D1_5a <- data.frame(FPR = rf_fpr_gb_0_5, TPR = rf_tpr_gb_0_5)
D1_5a$ID <- "RF_GB_1"
D2_5a <- data.frame(FPR = rf_fpr_gb_1_5, TPR = rf_tpr_gb_1_5)
D2_5a$ID <- "RF_GB_2"


D5_5 <- data.frame(FPR = fpr_full_0_5, TPR = tpr_full_0_5)
D5_5$ID <- "Full_1"
D8_5<- data.frame(FPR = fpr_full_4_5, TPR = tpr_full_4_5)
D8_5$ID <- "Full_2"

D5_5a <- data.frame(FPR = fpr_gb_0_5, TPR = tpr_gb_0_5)
D5_5a$ID <- "Full_GB_1"
D8_5a<- data.frame(FPR = fpr_gb_4_5, TPR = tpr_gb_4_5)
D8_5a$ID <- "Full_GB_2"

#Merge data frame

AUC5 <- c(
  trapz(D5_5$FPR,D5_5$TPR),
  trapz(D8_5$FPR,D8_5$TPR),
  trapz(D1_5$FPR,D1_5$TPR),
  trapz(D2_5$FPR,D2_5$TPR),
  trapz(D5_5a$FPR,D5_5a$TPR),
  trapz(D8_5a$FPR,D8_5a$TPR),
  trapz(D1_5a$FPR,D1_5a$TPR),
  trapz(D2_5a$FPR,D2_5a$TPR))

AUC5 <- as.data.frame(AUC5)
rownames(AUC5) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")

AUC1$ID <- factor(rownames(AUC1), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
AUC2$ID <- factor(rownames(AUC2), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
AUC3$ID <- factor(rownames(AUC3), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
AUC4$ID <- factor(rownames(AUC4), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
AUC5$ID <- factor(rownames(AUC5), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )


#Now do line plots
p1<-ggplot(data=AUC1, aes(x=ID,y=AUC1,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1)) 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_1.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p1)

p2<-ggplot(data=AUC2, aes(x=ID,y=AUC2,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_2.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p2)

p3<-ggplot(data=AUC3, aes(x=ID,y=AUC3,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_3.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p3)

p4<-ggplot(data=AUC4, aes(x=ID,y=AUC4,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_4.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p4)

p5<-ggplot(data=AUC5, aes(x=ID,y=AUC5,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_5.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p5)

saveRDS(AUC1,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_AUC_GB_1.rds",sep=""))
saveRDS(AUC2,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_AUC_GB_2.rds",sep=""))
saveRDS(AUC3,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_AUC_GB_3.rds",sep=""))
saveRDS(AUC4,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_AUC_GB_4.rds",sep=""))
saveRDS(AUC5,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_AUC_GB_5.rds",sep=""))



#Do plotting
d1 <- rbind(D5_1,D8_1,D1_1,D2_1,D5_1a,D8_1a,D1_1a,D2_1a)
p1 <- ggplot(d1, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_1_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p1)

d2 <- rbind(D5_2,D8_2,D1_2,D2_2,D5_2a,D8_2a,D1_2a,D2_2a)
p2 <- ggplot(d2, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_2_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p2)

d3 <- rbind(D5_3,D8_3,D1_3,D2_3,D5_3a,D8_3a,D1_3a,D2_3a)
p3 <- ggplot(d3, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_3_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p3)

d4 <- rbind(D5_4,D8_4,D1_5,D2_4,D5_4a,D8_4a,D1_4a,D2_4a)
p4 <- ggplot(d4, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_4_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p4)

d5 <- rbind(D5_5,D8_5,D1_5,D2_5,D5_5a,D8_5a,D1_5a,D2_5a)
p5 <- ggplot(d5, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC_GB_5_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p5)

colnames(AUC1)[1] <- "AUC"
colnames(AUC2)[1] <- "AUC"
colnames(AUC3)[1] <- "AUC"
colnames(AUC4)[1] <- "AUC"
colnames(AUC5)[1] <- "AUC"
AUCav <- data.frame(AUC= (AUC1$AUC+AUC2$AUC+AUC3$AUC+AUC4$AUC+AUC5$AUC)/5)
rownames(AUCav) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")
AUCav$ID <- factor(rownames(AUCav), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )

#rownames(AUCav) <- rownames(AUC1)
AUCav$ID <- rownames(AUCav)

p6<-ggplot(data=AUCav, aes(x=ID,y=AUC,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() + coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUCaverage_endo.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p6)

#AUC1 <- rbind(D1_1,D2_1,D3_1,D4_1,D5_1,D6_1,D7_1,D8_1,D9_1,D10_1)
#AUC2 <- rbind(D1_2,D2_2,D3_2,D4_2,D5_2,D6_2,D7_2,D8_2,D9_2,D10_2)
#AUC3 <- rbind(D1_3,D2_3,D3_3,D4_3,D5_3,D6_3,D7_3,D8_3,D9_3,D10_3)
#AUC4 <- rbind(D1_4,D2_4,D3_4,D4_4,D5_4,D6_4,D7_4,D8_4,D9_4,D10_4)
#AUC5 <- rbind(D1_5,D2_5,D3_5,D4_5,D5_5,D6_5,D7_5,D8_5,D9_5,D10_5)

#AUCb <- rbind(AUC1,AUC2,AUC3,AUC4,AUC5)
#colnames(AUCb) <- c("RF1","RF2","Targ1","Targ2","Full1","Chr1","Exp1","Full2","Chr2","Exp2","SE_Full_1","SE_Chrom_1","SE_Exp_1","SE_Full_2","SE_Chrom_2","SE_Exp_2","SE_nTL_1","SE_nTL_2","SE_TL_1","SE_TL_2")
#write.csv(as.data.frame(AUCb), file=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/AUC.csv",sep=""))


#AUC1 <- rbind(D1_1,D2_1,D6_1,D9_1)
#AUC2 <- rbind(D1_2,D2_2,D6_2,D9_2)
#AUC3 <- rbind(D1_3,D2_3,D6_3,D9_3)
#AUC4 <- rbind(D1_4,D2_4,D6_4,D9_4)
#AUC5 <- rbind(D1_5,D2_5,D6_5,D9_5)


#Do plotting
#p1 <- ggplot(AUC1, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p2 <- ggplot(AUC2, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p3 <- ggplot(AUC3, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p4 <- ggplot(AUC4, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p5 <- ggplot(AUC5, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#h <- grid.arrange(p1, p2, p3, p4, p5, ncol=5)

#ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/ROC_2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE, h)


#AUC1 <- rbind(D11_1,D12_1,D14_1,D15_1,D17_1,D18_1,D19_1,D20_1)
#AUC2 <- rbind(D11_2,D12_2,D14_2,D15_2,D17_2,D18_2,D19_2,D20_2)
#AUC3 <- rbind(D11_3,D12_3,D14_3,D15_3,D17_3,D18_3,D19_3,D20_3)
#AUC4 <- rbind(D11_4,D12_4,D14_4,D15_4,D17_4,D18_4,D19_4,D20_4)
#AUC5 <- rbind(D11_5,D12_5,D14_5,D15_5,D17_5,D18_5,D19_5,D20_5)


#Do plotting
#library(ggplot2)
#p1 <- ggplot(AUC1, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p2 <- ggplot(AUC2, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p3 <- ggplot(AUC3, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p4 <- ggplot(AUC4, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p5 <- ggplot(AUC5, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#h <- grid.arrange(p1, p2, p3, p4, p5, ncol=5)

#ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/ROC_4.pdf",sep=""),width = 50, height = 10,limitsize = FALSE, h)


#AUC1 <- rbind(D3_1,D4_1,D5_1,D6_1,D7_1,D8_1,D9_1,D10_1)
#AUC2 <- rbind(D3_2,D4_2,D5_2,D6_2,D7_2,D8_2,D9_2,D10_2)
#AUC3 <- rbind(D3_3,D4_3,D5_3,D6_3,D7_3,D8_3,D9_3,D10_3)
#AUC4 <- rbind(D3_4,D4_4,D5_4,D6_4,D7_4,D8_4,D9_4,D10_4)
#AUC5 <- rbind(D3_5,D4_5,D5_5,D6_5,D7_5,D8_5,D9_5,D10_5)


#Do plotting
#library(ggplot2)
#p1 <- ggplot(AUC1, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p2 <- ggplot(AUC2, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p3 <- ggplot(AUC3, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p4 <- ggplot(AUC4, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#p5 <- ggplot(AUC5, aes(x=FPR, y=TPR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
#  scale_color_manual(values=c("#006d2c","#00441b","#045a8d","#023858","#b30000","#7f0000","#252525","#000000","#54278f","#3f007d"))

#h <- grid.arrange(p1, p2, p3, p4, p5, ncol=5)

#ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/ROC_4.pdf",sep=""),width = 50, height = 10,limitsize = FALSE, h)

#Now the bbase ...

rm(list = ls())
np <- import("numpy")


#Random forests
rf_pr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears0_1.npy")
rf_pr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears0_2.npy")
rf_pr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears0_3.npy")
rf_pr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears0_4.npy")
rf_pr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears0_5.npy")
rf_pr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears1_1.npy")
rf_pr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears1_2.npy")
rf_pr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears1_3.npy")
rf_pr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears1_4.npy")
rf_pr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_pr_kears1_5.npy")

rf_re_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears0_1.npy")
rf_re_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears0_2.npy")
rf_re_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears0_3.npy")
rf_re_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears0_4.npy")
rf_re_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears0_5.npy")
rf_re_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears1_1.npy")
rf_re_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears1_2.npy")
rf_re_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears1_3.npy")
rf_re_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears1_4.npy")
rf_re_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoChr/rf_re_kears1_5.npy")

#NNs
pr_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears0_1.npy")
pr_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears0_2.npy")
pr_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears0_3.npy")
pr_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears0_4.npy")
pr_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears0_5.npy")
pr_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears1_1.npy")
pr_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears1_2.npy")
pr_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears1_3.npy")
pr_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears1_4.npy")
pr_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears1_5.npy")
pr_full_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears2_1.npy")
pr_full_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears2_2.npy")
pr_full_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears2_3.npy")
pr_full_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears2_4.npy")
pr_full_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears2_5.npy")
pr_full_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears4_1.npy")
pr_full_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears4_2.npy")
pr_full_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears4_3.npy")
pr_full_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears4_4.npy")
pr_full_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears4_5.npy")
pr_full_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears5_1.npy")
pr_full_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears5_2.npy")
pr_full_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears5_3.npy")
pr_full_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears5_4.npy")
pr_full_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears5_5.npy")
pr_full_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears6_1.npy")
pr_full_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears6_2.npy")
pr_full_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears6_3.npy")
pr_full_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears6_4.npy")
pr_full_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/pr_kears6_5.npy")

re_full_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears0_1.npy")
re_full_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears0_2.npy")
re_full_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears0_3.npy")
re_full_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears0_4.npy")
re_full_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears0_5.npy")
re_full_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears1_1.npy")
re_full_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears1_2.npy")
re_full_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears1_3.npy")
re_full_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears1_4.npy")
re_full_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears1_5.npy")
re_full_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears2_1.npy")
re_full_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears2_2.npy")
re_full_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears2_3.npy")
re_full_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears2_4.npy")
re_full_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears2_5.npy")
re_full_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears4_1.npy")
re_full_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears4_2.npy")
re_full_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears4_3.npy")
re_full_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears4_4.npy")
re_full_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears4_5.npy")
re_full_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears5_1.npy")
re_full_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears5_2.npy")
re_full_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears5_3.npy")
re_full_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears5_4.npy")
re_full_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears5_5.npy")
re_full_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears6_1.npy")
re_full_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears6_2.npy")
re_full_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears6_3.npy")
re_full_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears6_4.npy")
re_full_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFull/re_kears6_5.npy")


#Random forests for GB
rf_pr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears0_1.npy")
rf_pr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears0_2.npy")
rf_pr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears0_3.npy")
rf_pr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears0_4.npy")
rf_pr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears0_5.npy")
rf_pr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears1_1.npy")
rf_pr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears1_2.npy")
rf_pr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears1_3.npy")
rf_pr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears1_4.npy")
rf_pr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_pr_kears1_5.npy")

rf_re_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears0_1.npy")
rf_re_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears0_2.npy")
rf_re_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears0_3.npy")
rf_re_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears0_4.npy")
rf_re_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears0_5.npy")
rf_re_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears1_1.npy")
rf_re_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears1_2.npy")
rf_re_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears1_3.npy")
rf_re_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears1_4.npy")
rf_re_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/rf_re_kears1_5.npy")


#NNs for GB
pr_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears0_1.npy")
pr_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears0_2.npy")
pr_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears0_3.npy")
pr_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears0_4.npy")
pr_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears0_5.npy")
pr_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears1_1.npy")
pr_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears1_2.npy")
pr_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears1_3.npy")
pr_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears1_4.npy")
pr_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears1_5.npy")
pr_gb_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears2_1.npy")
pr_gb_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears2_2.npy")
pr_gb_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears2_3.npy")
pr_gb_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears2_4.npy")
pr_gb_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears2_5.npy")
pr_gb_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears4_1.npy")
pr_gb_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears4_2.npy")
pr_gb_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears4_3.npy")
pr_gb_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears4_4.npy")
pr_gb_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears4_5.npy")
pr_gb_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears5_1.npy")
pr_gb_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears5_2.npy")
pr_gb_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears5_3.npy")
pr_gb_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears5_4.npy")
pr_gb_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears5_5.npy")
pr_gb_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears6_1.npy")
pr_gb_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears6_2.npy")
pr_gb_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears6_3.npy")
pr_gb_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears6_4.npy")
pr_gb_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/pr_kears6_5.npy")

re_gb_0_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears0_1.npy")
re_gb_0_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears0_2.npy")
re_gb_0_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears0_3.npy")
re_gb_0_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears0_4.npy")
re_gb_0_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears0_5.npy")
re_gb_1_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears1_1.npy")
re_gb_1_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears1_2.npy")
re_gb_1_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears1_3.npy")
re_gb_1_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears1_4.npy")
re_gb_1_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears1_5.npy")
re_gb_2_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears2_1.npy")
re_gb_2_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears2_2.npy")
re_gb_2_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears2_3.npy")
re_gb_2_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears2_4.npy")
re_gb_2_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears2_5.npy")
re_gb_4_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears4_1.npy")
re_gb_4_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears4_2.npy")
re_gb_4_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears4_3.npy")
re_gb_4_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears4_4.npy")
re_gb_4_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears4_5.npy")
re_gb_5_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears5_1.npy")
re_gb_5_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears5_2.npy")
re_gb_5_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears5_3.npy")
re_gb_5_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears5_4.npy")
re_gb_5_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears5_5.npy")
re_gb_6_1 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears6_1.npy")
re_gb_6_2 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears6_2.npy")
re_gb_6_3 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears6_3.npy")
re_gb_6_4 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears6_4.npy")
re_gb_6_5 <- np$load("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/EndoFullGB/re_kears6_5.npy")

#Merge data frame
D1_1 <- data.frame(PR = rf_pr_full_0_1, RE = rf_re_full_0_1)
D1_1$ID <- "RF_1"
D2_1 <- data.frame(PR = rf_pr_full_1_1, RE = rf_re_full_1_1)
D2_1$ID <- "RF_2"
D1_1a <- data.frame(PR = rf_pr_gb_0_1, RE = rf_re_gb_0_1)
D1_1a$ID <- "RF_GB_1"
D2_1a <- data.frame(PR = rf_pr_gb_1_1, RE = rf_re_gb_1_1)
D2_1a$ID <- "RF_GB_2"

#
D5_1 <- data.frame(PR = pr_full_0_1, RE = re_full_0_1)
D5_1$ID <- "Full_1"
D8_1 <- data.frame(PR = pr_full_4_1, RE = re_full_4_1)
D8_1$ID <- "Full_2"
D5_1a <- data.frame(PR = pr_gb_0_1, RE = re_gb_0_1)
D5_1a$ID <- "Full_GB_1"
D8_1a <- data.frame(PR = pr_gb_4_1, RE = re_gb_4_1)
D8_1a$ID <- "Full_GB_2"


PR1 <- c(trapz(1-D5_1$RE,D5_1$PR),
          trapz(1-D8_1$RE,D8_1$PR),
          trapz(1-D1_1$RE,D1_1$PR),
          trapz(1-D2_1$RE,D2_1$PR),
          trapz(1-D5_1a$RE,D5_1a$PR),
          trapz(1-D8_1a$RE,D8_1a$PR),
          trapz(1-D1_1a$RE,D1_1a$PR),
          trapz(1-D2_1a$RE,D2_1a$PR))

PR1 <- as.data.frame(PR1)
rownames(PR1) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                    "GB1","GB2","RF1 (GB)","RF2 (GB)")


D1_2 <- data.frame(PR = rf_pr_full_0_2, RE = rf_re_full_0_2)
D1_2$ID <- "RF_1"
D2_2 <- data.frame(PR = rf_pr_full_1_2, RE = rf_re_full_1_2)
D2_2$ID <- "RF_2"
D1_2a <- data.frame(PR = rf_pr_gb_0_2, RE = rf_re_gb_0_2)
D1_2a$ID <- "RF_GB_1"
D2_2a <- data.frame(PR = rf_pr_gb_1_2, RE = rf_re_gb_1_2)
D2_2a$ID <- "RF_GB_2"


D5_2 <- data.frame(PR = pr_full_0_2, RE = re_full_0_2)
D5_2$ID <- "Full_1"
D8_2 <- data.frame(PR = pr_full_4_2, RE = re_full_4_2)
D8_2$ID <- "Full_2"
D5_2a <- data.frame(PR = pr_gb_0_2, RE = re_gb_0_2)
D5_2a$ID <- "Full_GB_1"
D8_2a <- data.frame(PR = pr_gb_4_2, RE = re_gb_4_2)
D8_2a$ID <- "Full_GB_2"

PR2 <- c(trapz(1-D5_2$RE,D5_2$PR),
         trapz(1-D8_2$RE,D8_2$PR),
         trapz(1-D1_2$RE,D1_2$PR),
         trapz(1-D2_2$RE,D2_2$PR),
         trapz(1-D5_2a$RE,D5_2a$PR),
         trapz(1-D8_2a$RE,D8_2a$PR),
         trapz(1-D1_2a$RE,D1_2a$PR),
         trapz(1-D2_2a$RE,D2_2a$PR))

PR2 <- as.data.frame(PR2)
rownames(PR2) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                   "GB1","GB2","RF1 (GB)","RF2 (GB)")



D1_3 <- data.frame(PR = rf_pr_full_0_3, RE = rf_re_full_0_3)
D1_3$ID <- "RF_1"
D2_3 <- data.frame(PR = rf_pr_full_1_3, RE = rf_re_full_1_3)
D2_3$ID <- "RF_2"
D1_3a <- data.frame(PR = rf_pr_gb_0_3, RE = rf_re_gb_0_3)
D1_3a$ID <- "RF_GB_1"
D2_3a <- data.frame(PR = rf_pr_gb_1_3, RE = rf_re_gb_1_3)
D2_3a$ID <- "RF_GB_2"


D5_3 <- data.frame(PR = pr_full_0_3, RE = re_full_0_3)
D5_3$ID <- "Full_1"
D8_3 <- data.frame(PR = pr_full_4_3, RE = re_full_4_3)
D8_3$ID <- "Full_2"
D5_3a <- data.frame(PR = pr_gb_0_3, RE = re_gb_0_3)
D5_3a$ID <- "Full_GB_1"
D8_3a <- data.frame(PR = pr_gb_4_3, RE = re_gb_4_3)
D8_3a$ID <- "Full_GB_2"

PR3 <- c(trapz(1-D5_3$RE,D5_3$PR),
         trapz(1-D8_3$RE,D8_3$PR),
         trapz(1-D1_3$RE,D1_3$PR),
         trapz(1-D2_3$RE,D2_3$PR),
         trapz(1-D5_3a$RE,D5_3a$PR),
         trapz(1-D8_3a$RE,D8_3a$PR),
         trapz(1-D1_3a$RE,D1_3a$PR),
         trapz(1-D2_3a$RE,D2_3a$PR))

PR3 <- as.data.frame(PR3)
rownames(PR3) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                   "GB1","GB2","RF1 (GB)","RF2 (GB)")


D1_4 <- data.frame(PR = rf_pr_full_0_4, RE = rf_re_full_0_4)
D1_4$ID <- "RF_1"
D2_4 <- data.frame(PR = rf_pr_full_1_4, RE = rf_re_full_1_4)
D2_4$ID <- "RF_2"
D1_4a <- data.frame(PR = rf_pr_gb_0_4, RE = rf_re_gb_0_4)
D1_4a$ID <- "RF_GB_1"
D2_4a <- data.frame(PR = rf_pr_gb_1_4, RE = rf_re_gb_1_4)
D2_4a$ID <- "RF_GB_2"

D5_4 <- data.frame(PR = pr_full_0_4, RE = re_full_0_4)
D5_4$ID <- "Full_1"
D8_4 <- data.frame(PR = pr_full_4_4, RE = re_full_4_4)
D8_4$ID <- "Full_2"

D5_4a <- data.frame(PR = pr_gb_0_4, RE = re_gb_0_4)
D5_4a$ID <- "Full_GB_1"
D8_4a <- data.frame(PR = pr_gb_4_4, RE = re_gb_4_4)
D8_4a$ID <- "Full_GB_2"

PR4 <- c(trapz(1-D5_4$RE,D5_4$PR),
         trapz(1-D8_4$RE,D8_4$PR),
         trapz(1-D1_4$RE,D1_4$PR),
         trapz(1-D2_4$RE,D2_4$PR),
         trapz(1-D5_4a$RE,D5_4a$PR),
         trapz(1-D8_4a$RE,D8_4a$PR),
         trapz(1-D1_4a$RE,D1_4a$PR),
         trapz(1-D2_4a$RE,D2_4a$PR))

PR4 <- as.data.frame(PR4)
rownames(PR4) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                   "GB1","GB2","RF1 (GB)","RF2 (GB)")



D1_5 <- data.frame(PR = rf_pr_full_0_5, RE = rf_re_full_0_5)
D1_5$ID <- "RF_1"
D2_5 <- data.frame(PR = rf_pr_full_1_5, RE = rf_re_full_1_5)
D2_5$ID <- "RF_2"

D1_5a <- data.frame(PR = rf_pr_gb_0_5, RE = rf_re_gb_0_5)
D1_5a$ID <- "RF2_GB_1"
D2_5a <- data.frame(PR = rf_pr_gb_1_5, RE = rf_re_gb_1_5)
D2_5a$ID <- "RF2_GB_2"

D5_5 <- data.frame(PR = pr_full_0_5, RE = re_full_0_5)
D5_5$ID <- "Full_1"
D8_5<- data.frame(PR = pr_full_4_5, RE = re_full_4_5)
D8_5$ID <- "Full_2"


D5_5a <- data.frame(PR = pr_gb_0_5, RE = re_gb_0_5)
D5_5a$ID <- "Full_GB_1"
D8_5a<- data.frame(PR = pr_gb_4_5, RE = re_gb_4_5)
D8_5a$ID <- "Full_GB_2"

PR5 <- c(trapz(1-D5_5$RE,D5_5$PR),
         trapz(1-D8_5$RE,D8_5$PR),
         trapz(1-D1_5$RE,D1_5$PR),
         trapz(1-D2_5$RE,D2_5$PR),
         trapz(1-D5_5a$RE,D5_5a$PR),
         trapz(1-D8_5a$RE,D8_5a$PR),
         trapz(1-D1_5a$RE,D1_5a$PR),
         trapz(1-D2_5a$RE,D2_5a$PR))

PR5 <- as.data.frame(PR5)
rownames(PR5) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                   "GB1","GB2","RF1 (GB)","RF2 (GB)")



PR1$ID <- factor(rownames(PR1), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
PR2$ID <- factor(rownames(PR2), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
PR3$ID <- factor(rownames(PR3), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
PR4$ID <- factor(rownames(PR4), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
PR5$ID <- factor(rownames(PR5), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )


p1<-ggplot(data=PR1, aes(x=ID,y=PR1,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_1.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p1)

p2<-ggplot(data=PR2, aes(x=ID,y=PR2,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_2.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p2)

p3<-ggplot(data=PR3, aes(x=ID,y=PR3,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_3.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p3)

p4<-ggplot(data=PR4, aes(x=ID,y=PR4,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_4.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p4)

p5<-ggplot(data=PR5, aes(x=ID,y=PR5,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() 
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_5.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p5)

 

saveRDS(PR1,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_PR_GB_1.rds",sep=""))
saveRDS(PR2,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_PR_GB_2.rds",sep=""))
saveRDS(PR3,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_PR_GB_3.rds",sep=""))
saveRDS(PR4,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_PR_GB_4.rds",sep=""))
saveRDS(PR5,paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/Endo_PR_GB_5.rds",sep=""))


#Do plotting
d1 <- rbind(D5_1,D8_1,D1_1,D2_1,D5_1a,D8_1a,D1_1a,D2_1a)
p1 <- ggplot(d1, aes(x=RE, y=PR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_1_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p1)

d2 <- rbind(D5_2,D8_2,D1_2,D2_2,D5_2a,D8_2a,D1_2a,D2_2a)
p2 <- ggplot(d2, aes(x=RE, y=PR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_2_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p2)

d3 <- rbind(D5_3,D8_3,D1_3,D2_3,D5_3a,D8_3a,D1_3a,D2_3a)
p3 <- ggplot(d3, aes(x=RE, y=PR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_3_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p3)

d4 <- rbind(D5_4,D8_4,D1_4,D2_4,D5_4a,D8_4a,D1_4a,D2_4a)
p4 <- ggplot(d4, aes(x=RE, y=PR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
 scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_4_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p4)

d5 <- rbind(D5_5,D8_5,D1_5,D2_5,D5_5a,D8_5a,D1_5a,D2_5a)
p5 <- ggplot(d5, aes(x=RE, y=PR, group=ID, color=ID)) + geom_line(aes(linetype=ID, color=ID)) + theme_classic() + scale_linetype_manual(values=c("solid","solid","dashed","dashed","longdash","longdash","dotted","dotted","twodash","twodash")) + 
  scale_color_manual(values=c("#023858","#023858","#0570b0","#0570b0","#045a8d","#045a8d","#3690c0","#3690c0"))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PR_GB_5_LP.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p5)


colnames(PR1)[1] <- "PR"
colnames(PR2)[1] <- "PR"
colnames(PR3)[1] <- "PR"
colnames(PR4)[1] <- "PR"
colnames(PR5)[1] <- "PR"
PRav <- data.frame(PR= (PR1$PR+PR2$PR+PR3$PR+PR4$PR+PR5$PR)/5)
rownames(PRav) <- c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)",
                   "GB1","GB2","RF1 (GB)","RF2 (GB)")
PRav$ID <- factor(rownames(PRav), levels = c("Prom1","Prom2","RF1 (Prom)","RF2 (Prom)","GB1","GB2","RF1 (GB)","RF2 (GB)") )
#rownames(PRav) <- rownames(PR1)
PRav$ID <- rownames(PRav)

p6<-ggplot(data=PRav, aes(x=ID,y=PR,fill=ID)) + geom_bar(stat="identity") + scale_fill_manual(values=c("#023858","#023858","#045a8d","#045a8d","#0570b0","#0570b0","#3690c0","#3690c0","#74a9cf","#74a9cf","#a6bddb","#a6bddb","#d0d1e6","#d0d1e6","#ece7f2","#fff7fb")) + theme_minimal() #+ coord_cartesian(ylim=c(0.5,1))
ggsave(filename=paste("/Volumes/GoogleDrive/My\ Drive/DeepXen/DeepXenJamboree/ExtendedROC/PRaverage_endo.pdf",sep=""),width = 15, height = 10,limitsize = FALSE, p6)

