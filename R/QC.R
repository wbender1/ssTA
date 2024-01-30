#############################################################################################################################################################################################################################################################
##QC
#############################################################################################################################################################################################################################################################
#' Quality Control
#'
#' @param fasta_file a fasta file containing nucleic acid sequences
#' @param seq_length_vector a vector containing the lengths of nucleic acid sequences in the fasta file
#' @param output_name name for output fasta file
#'
#' @return a list containing cleaned sequences is stored in the RStudio environment &  a fasta file containing cleaned sequences is written to the working directory.
#' @export
#'
#' @description Performs Quality Control measures to nucleic acid sequence data.
#' @description Pulls sequences of given lengths. Removes any sequence which has a character or symbol other than "ATCG".
#' @description Required package: seqinr
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples myLengths <- c(1725)
#' @examples myQCSeqs <- QC("PracticeSeqs_NA_GOOD.fasta", myLengths, "Output_GOOD/PracticeSeqs_Cleaned_GOOD.fasta")
#' @examples myLengths <- c(1725, 1726, 1727, 1729)
#' @examples myQCSeqs <- QC("PracticeSeqs_NA_MultiLengths.fasta", myLengths, "Output_MultiLengths/PracticeSeqs_Cleaned_MultiLengths.fasta")
#'
QC<- function(fasta_file, seq_length_vector, output_name){
  library(seqinr)
  mySeqs<- read.alignment(file = fasta_file, format = "fasta", forceToLower = FALSE) ##Read in fasta
  fastafile<- read.fasta(file = fasta_file, forceDNAtolower = FALSE) ##Read in fasta
  table(sapply(mySeqs$seq, nchar)) ##Check seq length
  plot(table(sapply(mySeqs$seq, nchar))) ##Check seq length
  numChar<-sapply(mySeqs$seq, nchar) ##Record the length of each seq and store as numChar
  indices_to_pull <- seq_length_vector
  pullMeAll <- numeric(0)
  for (i in indices_to_pull) {
    pullMe <- grep(i, numChar)
    pullMeAll <- c(pullMeAll, pullMe)
  }
  QC1file<-fastafile[pullMeAll] ##subset sequences from fastafile
  #print(names(QC1file))
  characters_to_search <- c("A", "T", "C", "G")
  namList <- character()
  for (i in 1:length(QC1file)) {
    myName <- names(QC1file[i])
    aSeq <- QC1file[[i]]
    if (!all(strsplit(aSeq, "")[[1]] %in% characters_to_search)) {
      namList <- c(namList, myName)
    }
  }
  QC2file<-QC1file[c(which(!(names(QC1file) %in% namList)))]
  write.fasta(as.list(QC2file), names(QC2file), file= output_name) ##Writes a new fasta file containing the seqs with correct length of NA
  return(QC2file)
}
