#############################################################################################################################################################################################################################################################
##Trait Association
#############################################################################################################################################################################################################################################################
#' Trait Association
#'
#' @param myPheno phenotypic data frame that contains a trait to be analyzed
#' @param col_name the column name of the trait to be analyzed
#' @param file_path directory path for the stats.csv file
#' @param file_name name of the stats.csv file
#'
#' @return a .csv file containing results from statistical tests
#' @export
#'
#' @examples read.csv(myPheno.csv)
#' @examples TA(myPheno, Disease_Severity, myStats, /myDirectory/)
TA<- function(myPheno, col_name, file_path, file_name){
  ##myPheno is a data frame that contains the trait we are measuring association of
  ##file_name is the name of the output stats .csv
  ##file_path is the working directory of where PNG files should be written to
  library(Rtsne)
  library(vegan)
  myPheno$col_name<-as.character(myPheno$col_name) ##should be character vector, sets Name/ID col as character
  myStats<-as.data.frame(matrix(data = NA, nrow = length(dir()), ncol = 5)) ##creates empty matrix to store stats in
  colnames(myStats)<-c("anosim.pval", "adonis2.pval", "bDis.perm.pval", "bDis.anova.pval","Names") ##Sets column names for myStats

  for(i in 1:length(dir())){
    myDM<-read.csv(dir()[i], row.names = 1, header = TRUE) ##Loop for every .csv in directory
    myDM1<-as.matrix(myDM) ##place myDM into myDM1 for later
    set.seed(1)
    myDM$Names<-rownames(myDM)
    intersect(myPheno$Global.Unique.Id, myDM$Names)
    allData2<-merge(myPheno, myDM, by.x = "Global.Unique.Id", by.y="Names")
    #Anosim
    myFac<-allData2$Corrected.RSV.Severity ##pull the column from allData2 of the trait you are testing association with and make a character vector from it
    myAnosim<-anosim(myDM1,myFac) ##perform anosim analysis on myDM1 and myFac
    myAnosim.pval<-myAnosim$signif ##pull significance number
    png(filename = paste0(file_path, dir()[i], "_AnosimPlot.png"), width = 7, height = 7, units = "in", pointsize = 12, bg = "white", res = 1000, type = "cairo")
    plot(myAnosim)
    dev.off()
    #Adonis
    myAdonis<-adonis2(myDM1~myFac) ##ANOVA using matrices, adonis2 statistical test on myDM1 and myFac
    myAdonis.pval<-myAdonis$`Pr(>F)`[1]
    #BetaDipsersion
    myDist<-as.dist(myDM1, diag = FALSE, upper = FALSE)
    dispersion<-betadisper(myDist, group=myFac) ##Implements Marti Anderson's PERMDISP2 procedure for the analysis of multivariate homogeneity of group dispersions (variances)
    #bDis perm #ANOVA like permutation test for Constrained Correspondence Analysis (cca), Redundancy Analysis (rda) or distance-based Redundancy Analysis (dbRDA, capscale) to assess the significance of constraints.
    bDis.perm<-permutest(dispersion)
    bDis.perm$tab$`Pr(>F)`[1]
    bDis.perm.pval<-bDis.perm$tab$`Pr(>F)`[1]
    #bDis.avova
    bDis.anova<-anova(betadisper(myDist, myFac)) ##Compute analysis of variance (or deviance) tables for one or more fitted model objects.
    bDis.anova
    bDis.anova.pval<-bDis.anova$`Pr(>F)`[1]
    png(filename = paste0(file_path, dir()[i], "_Hallrsv_bDisper_sdEllipse.png"), width = 7, height = 7, units = "in", pointsize = 12, bg = "white", res = 1000, type = "cairo")
    try(plot(dispersion, hull=FALSE, ellipse=TRUE)) ##sd ellipse
    dev.off()
    myStats[i,1:4]<-c(myAnosim.pval, myAdonis.pval, bDis.perm.pval, bDis.anova.pval)
    myStats[i,5]<-dir()[i]
    print(dir()[i])
  }
  write.csv(myStats, paste0(file_path, "TA_Stats_", file_name, "_.csv"))
  return(myStats)
}
