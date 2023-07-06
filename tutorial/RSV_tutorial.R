#Tutorial with RSV Sequences

library("ssTA")

setwd("~/ssTA/tutorial/")
QC("ssTA_tutorialSeqs.fasta", 1725, "~/ssTA/tutorial/QC_Seqs.fasta")

Align("QC_Seqs.fasta", "Aligned_Seqs.fasta")

DM("Aligned_Seqs.fasta", 575, "~/ssTA/tutorial/DM/tutorial")

setwd("~/ssTA/tutorial/DM/")
SBM("tutorial_DM.csv", "~/ssTA/tutorial/mySBM")

myPheno<-read.csv("~/ssTA/tutorial/myPheno.csv")
TA(myPheno$subtype, "~/ssTA/tutorial/tutorial", "tutorial", "tutorial")
