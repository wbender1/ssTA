#############################################################################################################################################################################################################################################################
##Trait Association
#############################################################################################################################################################################################################################################################
#' Trait Association
#'
#' @param trait_col column name which contains a trait to be analyzed
#' @param output_path directory path for output files
#' @param Plot_name name for output plots
#' @param Stats_name name for output stats file
#'
#' @return A .csv stats file containing anosim, adonis2, bDis.perm, and bDis.anova pvalues. A pdf file containing an Anosim plot. A pdf file containing a bDisper sdEllipse.
#' @export
#'
#' @description Performs 4 statistical tests (anosim, adonis2, bDis.perm, and bDis.anova) for association between a trait and a distance matrix.
#' @description Required packages: Rtsne, vegan
#'
#' @examples setwd("~/ssTA/tutorial/")
#' @examples myPheno<-read.csv("myPheno.csv")
#' @examples setwd("DM_GOOD")
#' @examples myStats <- TA(2, "../Output_GOOD/", "PracticeSeqs_GOOD", "PracticeSeqs_GOOD")
#' @examples setwd("..")

TA<-function(trait_col_num, output_path, Plot_name, Stats_name){
myStats<-as.data.frame(matrix(data = NA, nrow = length(dir()), ncol = 5)) ##creates empty matrix to store stats in
colnames(myStats)<-c("anosim.pval", "adonis2.pval", "bDis.perm.pval", "bDis.anova.pval","Names") ##Sets column names for myStats
library(Rtsne)
library(vegan)
for(i in 1:length(dir())){
  myDM<-read.csv(dir()[i], row.names = 1)
  myDM1<-as.matrix(myDM) ##place myDM into myDM1 for later
  set.seed(1) # for reproducibility
  myDM$Names<-rownames(myDM) ##create a column "Names" in myDM and then adds the row names from myDM into the "Names" column
  myPhenoSubIndices <- which(myPheno$name %in% myDF$names)
  myPhenoSub <- myPheno[myPhenoSubIndices,]
  intersect(myDM$Names, myPhenoSub$name) ##check to see that the Global.Unique.Id column matches with the names found in myDM$Names
  allData2<-merge(myDM, myPhenoSub, by.x = "Names", by.y="name") ##merge the data by matching the Global.Unique.Id column and myDM$Names and push into allData2 data frame
#Anosim
  myFac<-myPhenoSub[,trait_col_num] ##pull the Corrected.RSV.Severity column from allData2 and make a character vector with it
  myAnosim<-anosim(myDM1,myFac) ##perform anosim analysis on myDM1 and myFac
  myAnosim.pval<-myAnosim$signif ##pull significance number and send to myAnosim.pval
  pdf(paste0(output_path, "AnosimPlot_", Plot_name,".pdf")) ##set the directory path for the directory which you want to write Anosim plots to and open a pdf
  plot(myAnosim) ##plot myAnosim to the pdf
  dev.off() ##close the pdf
#Adonis
  myAdonis<-adonis2(myDM1~myFac) ##ANOVA using matrices, adonis2 statistical test on myDM1 and myFac
  myAdonis.pval<-myAdonis$`Pr(>F)`[1] ##
#BetaDipsersion
#Non-euclidean distances between objects and group centroids are handled by reducing the original distances to principal coordinates.
  myDist<-as.dist(myDM1, diag = FALSE, upper = FALSE) ##
  dispersion<-betadisper(myDist, group=myFac) ##Implements Marti Anderson's PERMDISP2 procedure for the analysis of multivariate homogeneity of group dispersions (variances)
#bDis perm #ANOVA like permutation test for Constrained Correspondence Analysis (cca), Redundancy Analysis (rda) or distance-based Redundancy Analysis (dbRDA, capscale) to assess the significance of constraints.
  bDis.perm<-permutest(dispersion) ##
  bDis.perm$tab$`Pr(>F)`[1] ##
  bDis.perm.pval<-bDis.perm$tab$`Pr(>F)`[1] ##
#bDis.avova
  bDis.anova<-anova(betadisper(myDist, myFac)) ##Compute analysis of variance (or deviance) tables for one or more fitted model objects.
  bDis.anova
  bDis.anova.pval<-bDis.anova$`Pr(>F)`[1] ##
  pdf(paste0(output_path, "bDisper_sdEllipse_", Plot_name, ".pdf")) ##set the directory path for the directory which you want to write Ellipse plots to and open a pdf
  try(plot(dispersion, hull=FALSE, ellipse=TRUE)) ##sd ellipse
  dev.off() ##close pdf
  myStats[i,1:4]<-c(myAnosim.pval, myAdonis.pval, bDis.perm.pval, bDis.anova.pval) ##place the statistics produced in loop into myStats
  myStats[i,5]<-Plot_name ##place the name of each file into the 5th column of myStats
  print(Plot_name) ##print file name to directory
}
write.csv(myStats, paste0(output_path, "HDstats_", Stats_name,".csv"))
return(myStats)
}


