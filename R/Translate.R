#############################################################################################################################################################################################################################################################
##Translate
#############################################################################################################################################################################################################################################################
#' Translate
#'
#' @param input_file name of a fasta file containing Nucleic Acid sequences
#' @param output_file name of the output fasta file containing Amino Acid sequences
#'
#' @return a fasta file containing Amino acid sequences
#' @export
#'
#' @description Translates Nucleic acid sequences using bioseq seq_translate
#' @description Required package: bioseq
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples myTranslatedSeqs <- Translate("Output_GOOD/PracticeSeqs_Cleaned_GOOD.fasta", "Output_GOOD/PracticeSeqs_Translated_GOOD.fasta")
#'
Translate<- function(input_file, output_file){
  library(bioseq)
  mySeqs1<-read.alignment(input_file, format = "fasta", forceToLower = FALSE)
  timestamp1 <- timestamp(); time1 <- as.POSIXct(strptime(timestamp1, "##------ %a %b %e %H:%M:%S %Y ------##"))
  mydna<-lapply(mySeqs1$seq, dna)
  AA<-lapply(mydna, seq_translate)
  timestamp2 <- timestamp(); time2 <- as.POSIXct(strptime(timestamp2, "##------ %a %b %e %H:%M:%S %Y ------##"))
  write.fasta(AA, names = mySeqs1$nam, file= output_file)
  time_difference <- as.numeric(difftime(time2, time1, units = "secs")); cat("Time difference: ", time_difference/60, " minutes\n")
  return(AA)
  PRINT(time_difference)
}
