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
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples Align("CleanedSeqs.fasta", "myAlignedSeqs.fasta")
#'
Align<- function(input_file, output_file){
  library(msa)
  orf <- readAAStringSet(file = input_file, format = "fasta") ##set name of fasta file containing the seqs with correct length/ previously translated to AA
  timestamp()
  musAlign<-msa(orf, "ClustalOmega") ##perform msa to align sequences
  timestamp()
  writeXStringSet(unmasked(musAlign), file= output_file) ##set name of file
}
