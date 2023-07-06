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
#' @examples SBM(myMatrix, /myDirectory/myFileName)
SBM<-function(input_file, output_file){
  library(stringdist)
  ##Dimension reduction to allow for 2D plotting
  ##Classical multidimensional scaling on distance matrix
  myMDS<-cmdscale(input_file, k=2, eig=T)
  ##Check goodness of fit
  myMDS$GOF
  ##Plot the reduced data
  pdf(file = paste0(output_file, "SBM_plot.pdf"))
  plot(myMDS$points[,1],myMDS$points[,2])
  #text(myMDS$points[,1],myMDS$points[,2], labels=rownames(myMDS$points), cex = .5)
  dev.off()
  myDF<-as.data.frame(myMDS["points"])
  myDF$names<-rownames(myDF)
  write.csv(myDF, file = output_file)
}
