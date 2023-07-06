#############################################################################################################################################################################################################################################################
##Trait Association
#############################################################################################################################################################################################################################################################
#' Trait Association
#'
#' @param trait_col column name which contains a trait to be analyzed
#' @param path directory path for output files
#' @param Plot_name name for output plots
#' @param Stats_name name for output stats file
#'
#' @return A .csv stats file containing anosim, adonis2, bDis.perm, and bDis.anova pvalues. A pdf file containing an Anosim plot. A pdf file containing a bDisper sdEllipse.
#' @export
#'
#' @description Performs 4 statistical tests for association between a trait and a distance matrix.
#'
#' @examples setwd("~/ssTA/tutorial/DM/") ##set the directory path for the directory which contains all of your Distance Matrix files
#' @examples myPheno<-read.csv("~/ssTA/tutorial/myPheno.csv") ##Read in metadata file
#' @examples myPheno$name<- as.character(myPheno$name) ##should be character vector, sets GenBank.Accession col as character
#' @examples TA(myPheno$subtype, "~/ssTA/tutorial/", "Test", "Test")

TA<-function(trait_col, path, Plot_name, Stats_name){
myStats<-as.data.frame(matrix(data = NA, nrow = length(dir()), ncol = 5)) ##creates empty matrix to store stats in
colnames(myStats)<-c("anosim.pval", "adonis2.pval", "bDis.perm.pval", "bDis.anova.pval","Names") ##Sets column names for myStats

for(i in 1:length(dir())){
  myDM<-read.csv(dir()[i], row.names = 1)
  myDM1<-as.matrix(myDM) ##place myDM into myDM1 for later
  set.seed(1) # for reproducibility
  myDM$Names<-rownames(myDM) ##create a column "Names" in myDM and then adds the row names from myDM into the "Names" column
  intersect(myPheno$name, myDM$Names) ##check to see that the Global.Unique.Id column matches with the names found in myDM$Names
  allData2<-merge(myPheno, myDM, by.x = "name", by.y="Names") ##merge the data by matching the Global.Unique.Id column and myDM$Names and push into allData2 data frame
#Anosim
  myFac<-trait_col ##pull the Corrected.RSV.Severity column from allData2 and make a character vector with it
  myAnosim<-anosim(myDM1,myFac) ##perform anosim analysis on myDM1 and myFac
  myAnosim.pval<-myAnosim$signif ##pull significance number and send to myAnosim.pval
  pdf(paste0(path, Plot_name,"_AnosimPlot.pdf")) ##set the directory path for the directory which you want to write Anosim plots to and open a pdf
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
  pdf(paste0(path, Plot_name, "_bDisper_sdEllipse.pdf")) ##set the directory path for the directory which you want to write Ellipse plots to and open a pdf
  try(plot(dispersion, hull=FALSE, ellipse=TRUE)) ##sd ellipse
  dev.off() ##close pdf
  myStats[i,1:4]<-c(myAnosim.pval, myAdonis.pval, bDis.perm.pval, bDis.anova.pval) ##place the statistics produced in loop into myStats
  myStats[i,5]<-Plot_name ##place the name of each file into the 5th column of myStats
  print(Plot_name) ##print file name to directory
}

return(myStats)
write.csv(myStats, paste0(path, Stats_name,"_HDstats.csv"))
}


