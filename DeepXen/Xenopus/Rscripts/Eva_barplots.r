library(reticulate)
library(pheatmap)

#Load the endoderm counts
counts <- read.table("/Volumes/Overflow/Scores/EndoFull/Counts.csv",sep = ',', header = F)
counts <- counts[2:dim(counts)[1],2:dim(counts)[2]]
colnames(counts) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts) <- c("Full_seed1","Neg_seed1","Pos_seed1","Full_seed2","Neg_seed2","Pos_seed2")

#Now do TL?
counts2 <- read.table("/Volumes/Overflow/Scores/EndoFull/Count_TLs.csv",sep = ',', header = F)
counts2 <- counts2[2:3,2:dim(counts2)[2]]
colnames(counts2) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts2) <- c("TL_seed1","TL_seed2")

#Now do target only 
counts3 <- read.table("/Volumes/Overflow/Scores/EndoFull/Counts_endoTargOnly.csv",sep = ',', header = F)
counts3 <- counts3[2,2:dim(counts3)[2]]
colnames(counts3) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts3) <- c("Targ_seed1")

#Summarise the data
Endo <- rbind(counts,counts2,counts3)



#Now load SomEndo

#Now load Endo_GB

#ow load SomEndo_GB


#Load the endoderm counts
counts <- read.table("/Volumes/Overflow/Scores/SomEndoFull/Counts.csv",sep = ',', header = F)
counts <- counts[2:dim(counts)[1],2:dim(counts)[2]]
colnames(counts) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts) <- c("Full_seed1","Neg_seed1","Pos_seed1","Full_seed2","Neg_seed2","Pos_seed2")

#Now do TL?
counts2 <- read.table("/Volumes/Overflow/Scores/SomEndoFull/Count_TLs.csv",sep = ',', header = F)
counts2 <- counts2[2:3,2:dim(counts2)[2]]
colnames(counts2) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts2) <- c("TL_seed1","TL_seed2")

#Now do target only 
counts3 <- read.table("/Volumes/Overflow/Scores/SomEndoFull/Counts_ectoTargOnly.csv",sep = ',', header = F)
counts3 <- counts3[2,2:dim(counts3)[2]]
colnames(counts3) <- c("OnP","OnT","OnTP","OnFP","OnTN","OnFN","OffP","OffT","OffTP","OffFP","OnTN","OffFN","RUP","RUT","RUTP","RUFP","RUTN","RUFN","RDP","RDT","RDTP","RDFP","RDTN","RDFN","OtherP","OtherT")
rownames(counts3) <- c("Targ_seed1")

#Summarise the data
SomEcto <- rbind(counts,counts2,counts3)



