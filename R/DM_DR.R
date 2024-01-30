#############################################################################################################################################################################################################################################################
##Distance Matrix Dimension Reduction
#############################################################################################################################################################################################################################################################
#' Distance Matrix Dimension Reduction
#'
#' @param input_file a distance matrix
#' @param output_file the name and/or directory path for the output .pdf and .csv
#'
#' @return a .csv file containing a distance matrix and a .pdf file which contains the plotted points
#' @export
#'
#' @description Performs multiple dimensional scaling with CMDscale, writes a .csv file containing coordinates, and produces a plot that represents the antigenic distances between sequences.
#' @description Required package: stringdist
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples myDF <- DM_DR("DM_GOOD/PracticeSeqs_GOOD_DM.csv", "Output_GOOD/PracticeSeqs_points_GOOD")
#'
DM_DR<-function(input_file, output_file){
  library(stringdist)
  ##Dimension reduction to allow for 2D plotting
  ##Classical multidimensional scaling on distance matrix
  myMat<- as.matrix(read.csv(input_file, row.names = 1))
  myMDS<-cmdscale(myMat, k=2, eig=T)
  ##Check goodness of fit
  myMDS$GOF
  ##Plot the reduced data
  pdf(file = paste0(output_file, ".pdf"))
  plot(myMDS$points[,1],myMDS$points[,2])
  dev.off()
  myPlot <- plot(myMDS$points[,1],myMDS$points[,2])
  myDF<-as.data.frame(myMDS["points"])
  myDF$names<-rownames(myDF)
  write.csv(myDF, file= paste0(output_file, ".csv"))
  return(myDF)
  return(myPlot)
}
