#############################################################################################################################################################################################################################################################
##Alignment
#############################################################################################################################################################################################################################################################
#' Alignment
#'
#' @param input_file name of a fasta file containing Amino Acid sequences
#' @param output_file name of the output fasta file
#'
#' @return a fasta file containing aligned sequences
#' @export
#'
#' @description Aligns amino acid sequences using msa ClustalOmega.
#' @description Required package: msa
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples myAlignedSeqs <- Align("Output_GOOD/PracticeSeqs_Translated_GOOD.fasta", "Output_GOOD/PracticeSeqs_Aligned_GOOD.fasta")
#'
Align<- function(input_file, output_file){
  library(msa)
  orf <- readAAStringSet(file = input_file, format = "fasta") ##set name of fasta file containing the seqs with correct length/ previously translated to AA
  timestamp1 <- timestamp(); time1 <- as.POSIXct(strptime(timestamp1, "##------ %a %b %e %H:%M:%S %Y ------##"))
  musAlign<-msa(orf, "ClustalOmega") ##perform msa to align sequences
  timestamp2 <- timestamp(); time2 <- as.POSIXct(strptime(timestamp2, "##------ %a %b %e %H:%M:%S %Y ------##"))
  writeXStringSet(unmasked(musAlign), file= output_file) ##set name of file
  time_difference <- as.numeric(difftime(time2, time1, units = "secs")); cat("Time difference: ", time_difference/60, " minutes\n")
  return(musAlign)
  PRINT(time_difference)
}
