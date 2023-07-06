#############################################################################################################################################################################################################################################################
##Sequnce Based Mapping
#############################################################################################################################################################################################################################################################
#' Sequence Based Mapping
#'
#' @param input_file a distance matrix
#' @param output_file the name and/or directory path for the output .pdf and .csv
#'
#' @return a .csv file containing a distance matrix and a .pdf file which contains the plotted points
#' @export
#'
#' @description Performs multiple dimensional scaling with CMDscale and plots the antigenic distances between sequences.
#'
#' @examples setwd("~/ssTA/tutorial/DM/")
#' @examples SBM("tutorial_DM.csv", "~/ssTA/tutorial/tutorial_DM")
#'
SBM<-function(input_file, output_file){
  library(stringdist)
  ##Dimension reduction to allow for 2D plotting
  ##Classical multidimensional scaling on distance matrix
  myMat<- as.matrix(read.csv(input_file, row.names = 1))
  myMDS<-cmdscale(myMat, k=2, eig=T)
  ##Check goodness of fit
  myMDS$GOF
  ##Plot the reduced data
  pdf(file = paste0(output_file, "_plot.pdf"))
  plot(myMDS$points[,1],myMDS$points[,2])
  dev.off()
  myDF<-as.data.frame(myMDS["points"])
  myDF$names<-rownames(myDF)
  write.csv(myDF, file= paste0(output_file, "_points.csv"))
}
