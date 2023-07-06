#############################################################################################################################################################################################################################################################
##QC
#############################################################################################################################################################################################################################################################
#' Quality Control
#'
#' @param fasta_file a fasta file containing nucleic acid sequences
#' @param seq_length length of sequences
#' @param output_name name for output fasta file
#'
#' @return a fasta file containing cleaned sequences is written
#' @export
#'
#' @description Performs Quality Control measures to nucleic acid sequence data. Removes any sequence which has a letter other than "ATCG".
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples QC("Seqs.fasta", 1725, "CleanedSeqs.fasta")
#'
QC<- function(fasta_file, seq_length, output_name){
  library(seqinr)

  mySeqs<- read.alignment(file = fasta_file, format = "fasta", forceToLower = FALSE) ##Read in fasta
  fastafile<- read.fasta(file = fasta_file) ##Read in fasta
  table(sapply(mySeqs$seq, nchar)) ##Check seq length
  plot(table(sapply(mySeqs$seq, nchar))) ##Check seq length
  numChar<-sapply(mySeqs$seq, nchar) ##Record the length of each seq and store as numChar
  pullMeAll<-grep(seq_length, numChar) ##Pull the indices of every seq length that is equal to 195
  QC1file<-fastafile[pullMeAll] ##subset sequences from fastafile
  namList<-list()
  for (i in 1:length(QC1file)) {
    myName<-names(QC1file[i])
    aSeq<-QC1file[[i]]
    hasB<-which(strsplit(aSeq, "")[[1]]=="b")
    hasD<-which(strsplit(aSeq, "")[[1]]=="d")
    hasE<-which(strsplit(aSeq, "")[[1]]=="e")
    hasF<-which(strsplit(aSeq, "")[[1]]=="f")
    hasH<-which(strsplit(aSeq, "")[[1]]=="h")
    hasI<-which(strsplit(aSeq, "")[[1]]=="i")
    hasJ<-which(strsplit(aSeq, "")[[1]]=="j")
    hasK<-which(strsplit(aSeq, "")[[1]]=="k")
    hasL<-which(strsplit(aSeq, "")[[1]]=="l")
    hasM<-which(strsplit(aSeq, "")[[1]]=="m")
    hasN<-which(strsplit(aSeq, "")[[1]]=="n")
    hasO<-which(strsplit(aSeq, "")[[1]]=="o")
    hasP<-which(strsplit(aSeq, "")[[1]]=="p")
    hasQ<-which(strsplit(aSeq, "")[[1]]=="q")
    hasR<-which(strsplit(aSeq, "")[[1]]=="r")
    hasS<-which(strsplit(aSeq, "")[[1]]=="s")
    hasU<-which(strsplit(aSeq, "")[[1]]=="u")
    hasV<-which(strsplit(aSeq, "")[[1]]=="v")
    hasW<-which(strsplit(aSeq, "")[[1]]=="w")
    hasX<-which(strsplit(aSeq, "")[[1]]=="x")
    hasY<-which(strsplit(aSeq, "")[[1]]=="y")
    hasZ<-which(strsplit(aSeq, "")[[1]]=="z")
    if(length(hasB)>0){
      namList<-c(namList,myName)}
    if(length(hasD)>0){
      namList<-c(namList,myName)}
    if(length(hasE)>0){
      namList<-c(namList,myName)}
    if(length(hasF)>0){
      namList<-c(namList,myName)}
    if(length(hasH)>0){
      namList<-c(namList,myName)}
    if(length(hasI)>0){
      namList<-c(namList,myName)}
    if(length(hasJ)>0){
      namList<-c(namList,myName)}
    if(length(hasK)>0){
      namList<-c(namList,myName)}
    if(length(hasL)>0){
      namList<-c(namList,myName)}
    if(length(hasM)>0){
      namList<-c(namList,myName)}
    if(length(hasN)>0){
      namList<-c(namList,myName)}
    if(length(hasO)>0){
      namList<-c(namList,myName)}
    if(length(hasP)>0){
      namList<-c(namList,myName)}
    if(length(hasQ)>0){
      namList<-c(namList,myName)}
    if(length(hasR)>0){
      namList<-c(namList,myName)}
    if(length(hasS)>0){
      namList<-c(namList,myName)}
    if(length(hasU)>0){
      namList<-c(namList,myName)}
    if(length(hasV)>0){
      namList<-c(namList,myName)}
    if(length(hasW)>0){
      namList<-c(namList,myName)}
    if(length(hasX)>0){
      namList<-c(namList,myName)}
    if(length(hasY)>0){
      namList<-c(namList,myName)}
    if(length(hasZ)>0){
      namList<-c(namList,myName)}
  }
  QC2file<-QC1file[c(which(!(names(QC1file) %in% namList)))]
  QC3file<-QC1file[c(which(names(QC1file) %in% namList))]
  write.fasta(as.list(QC2file), names(QC2file), file= output_name) ##Writes a new fasta file containing the seqs with correct length of NA
}
