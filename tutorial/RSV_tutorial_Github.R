#ssTA Tutorials
library("ssTA")
setwd("~/ssTA_2.0/tutorial")
#############################################################################################################################################################################################################################################################
#Tutorial with GOOD Nucleic Acid sequences. All sequences pass through QC function.
#############################################################################################################################################################################################################################################################
#Create a vector containing all the lengths of sequences that are contained in the input fasta file.
myLengths <- c(1725)

#Perform QC to pull all sequences of previously recorded lengths and remove any sequences which has a character or symbol other than "ATCG".
myQCSeqs <- QC("PracticeSeqs_NA_GOOD.fasta", myLengths, "Output_GOOD/PracticeSeqs_Cleaned_GOOD.fasta")

#Perform Translate to translate Nucleic Acid sequences into Amino Acid sequences if necessary.
myTranslatedSeqs <- Translate("Output_GOOD/PracticeSeqs_Cleaned_GOOD.fasta", "Output_GOOD/PracticeSeqs_Translated_GOOD.fasta")

#Perform Align to align sequences to the same length if necessary.
myAlignedSeqs <- Align("Output_GOOD/PracticeSeqs_Translated_GOOD.fasta", "Output_GOOD/PracticeSeqs_Aligned_GOOD.fasta")

#Perform Distance Matrix to calculate the Antigenic Hamming Distances between our sequences and store the values as a distance matrix.
myDM <- DM("Output_GOOD/PracticeSeqs_Aligned_GOOD.fasta", 575, "DM_GOOD/PracticeSeqs_GOOD")

#Perform Distance Matrix Dimension Reduction to reduce the Distance Matrix into a set of coordinates.
myDF <- DM_DR("DM_GOOD/PracticeSeqs_GOOD_DM.csv", "Output_GOOD/PracticeSeqs_points_GOOD")

#Read in Phenotype data.
myPheno<-read.csv("myPheno.csv")

#Change to the directory which contains Distance Matrix/Matrices.
setwd("DM_GOOD")

#Perform Trait Association to calculate anosim, adonis2, bDis.perm, and bDis.anova pvalues.
myStats <- TA(2, "../Output_GOOD/", "PracticeSeqs_GOOD", "PracticeSeqs_GOOD")

#Return to parent directory
setwd("..")
#############################################################################################################################################################################################################################################################
#Tutorial with BAD Nucleic Acid sequences. Only 6 of the sequences pass through QC function.
#############################################################################################################################################################################################################################################################
#Create a vector containing all the lengths of sequences that are contained in the input fasta file.
myLengths <- c(1725)

#Perform QC to pull all sequences of previously recorded lengths and remove any sequences which has a character or symbol other than "ATCG".
myQCSeqs <- QC("PracticeSeqs_NA_BAD.fasta", myLengths, "Output_BAD/PracticeSeqs_Cleaned_BAD.fasta")

#Perform Translate to translate Nucleic Acid sequences into Amino Acid sequences if necessary.
myTranslatedSeqs <- Translate("Output_BAD/PracticeSeqs_Cleaned_BAD.fasta", "Output_BAD/PracticeSeqs_Translated_BAD.fasta")

#Perform Align to align sequences to the same length if necessary.
myAlignedSeqs <- Align("Output_BAD/PracticeSeqs_Translated_BAD.fasta", "Output_BAD/PracticeSeqs_Aligned_BAD.fasta")

#Perform Distance Matrix to calculate the Antigenic Hamming Distances between our sequences and store the values as a distance matrix.
myDM <- DM("Output_BAD/PracticeSeqs_Aligned_BAD.fasta", 575, "DM_BAD/PracticeSeqs_BAD")

#Perform Distance Matrix Dimension Reduction to reduce the Distance Matrix into a set of coordinates.
myDF <- DM_DR("DM_BAD/PracticeSeqs_BAD_DM.csv", "Output_BAD/PracticeSeqs_points_BAD")

#Read in Phenotype data.
myPheno<-read.csv("myPheno.csv")

#Change to the directory which contains Distance Matrix/Matrices.
setwd("DM_BAD")

#Perform Trait Association to calculate anosim, adonis2, bDis.perm, and bDis.anova pvalues.
myStats <- TA(2, "../Output_BAD/", "PracticeSeqs_BAD", "PracticeSeqs_BAD")

#Return to parent directory
setwd("..")
#############################################################################################################################################################################################################################################################
#Tutorial with GOOD Nucleic Acid sequences of multiple lengths. This will require pulling multiple length of sequences and aligning them to the same length.
#For the sake of this example, we want any sequence that has length greater than or equal to 1725 Nucleic Acids. 9 of the 10 sequences will be selected.
#############################################################################################################################################################################################################################################################
#Create a vector containing all the lengths of sequences that are contained in the input fasta file.
myLengths <- c(1725, 1726, 1727, 1729)

#Perform QC to pull all sequences of previously recorded lengths and remove any sequences which has a character or symbol other than "ATCG".
myQCSeqs <- QC("PracticeSeqs_NA_MultiLengths.fasta", myLengths, "Output_MultiLengths/PracticeSeqs_Cleaned_MultiLengths.fasta")

#Perform Translate to translate Nucleic Acid sequences into Amino Acid sequences if necessary.
myTranslatedSeqs <- Translate("Output_MultiLengths/PracticeSeqs_Cleaned_MultiLengths.fasta", "Output_MultiLengths/PracticeSeqs_Translated_MultiLengths.fasta")

#Perform Align to align sequences to the same length if necessary.
myAlignedSeqs <- Align("Output_MultiLengths/PracticeSeqs_Translated_MultiLengths.fasta", "Output_MultiLengths/PracticeSeqs_Aligned_MultiLengths.fasta")

#Perform Distance Matrix to calculate the Antigenic Hamming Distances between our sequences and store the values as a distance matrix.
myDM <- DM("Output_MultiLengths/PracticeSeqs_Aligned_MultiLengths.fasta", 575, "DM_MultiLengths/PracticeSeqs_MultiLengths")

#Perform Distance Matrix Dimension Reduction to reduce the Distance Matrix into a set of coordinates.
myDF <- DM_DR("DM_MultiLengths/PracticeSeqs_MultiLengths_DM.csv", "Output_MultiLengths/PracticeSeqs_points_MultiLengths")

#Read in Phenotype data.
myPheno<-read.csv("myPheno.csv")

#Change to the directory which contains Distance Matrix/Matrices.
setwd("DM_MultiLengths")

#Perform Trait Association to calculate anosim, adonis2, bDis.perm, and bDis.anova pvalues.
myStats <- TA(2, "../Output_MultiLengths/", "PracticeSeqs_MultiLengths", "PracticeSeqs_MultiLengths")

#Return to parent directory
setwd("..")

