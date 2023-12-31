#' Distance Matrix
#'
#' @param input_file a fasta file containing aligned amino acid sequences
#' @param sequence_length length of sequences in fasta file
#' @param output_file the name and/or directory path for the output .csv
#'
#' @return a .csv file containing a distance matrix of our sequences
#' @export
#'
#' @description Creates a distance matrix for a set of aligned sequences.
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples DM("myAlignedSeqs.fasta", 575, "~/ssTA/tutorial/DM/TEST")
#'
DM<-function(input_file, sequence_length, output_file){
  library(stringdist)
  viruAln<-read.alignment(file = input_file, format = "fasta", forceToLower = FALSE)
  table(sapply(viruAln$seq, nchar))
  plot(table(sapply(viruAln$seq, nchar)))
  for (i in 1:length(viruAln$nam)) {
    name<-viruAln$nam[[i]]
    seq<-viruAln$seq[[i]]
    assign(name, seq)
    #A epitope
    a<-substr(seq, 1, sequence_length)
    A<-paste0(a)
    assign(name, A)
    if(i == 1){
      namList<-name}
    else {namList<-c(namList,name)}
  }
  myMat<-matrix(data=NA, nrow=length(namList), ncol=length(namList))
  colnames(myMat)<-namList
  rownames(myMat)<-namList
  for (i in 1:length(namList)){
    name<-paste0(namList[i],"_vec")
    v1<-get(namList[i])
    for (j in 1:length(namList)){
      v2<-get(namList[j])
      dis<-stringdist(v1, v2, method="h")
      if(j==1){vec<-dis} else{vec<-c(vec, dis)}
      assign(name,vec)}
    if(j==length(namList)){myMat[, i]<-vec}
  }
  ##This results in a n x n matrix which also shows length of aligned sequences
  write.csv(myMat, file= paste0(output_file, "_DM.csv"))
}
