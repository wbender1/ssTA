
#############################################################################################################################################################################################################################################################
##QC
#############################################################################################################################################################################################################################################################
setwd("/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/seqData_raw")
QC<- function(X, Y, V, Z){
  library(seqinr)

  mySeqs<- read.alignment(file = X, format = "fasta", forceToLower = FALSE) ##Read in fasta
  fastafile<- read.fasta(file = X) ##Read in fasta
  myNames<-gsub('-', '', mySeqs$nam)
  table(sapply(mySeqs$seq, nchar)) ##Check seq length
  plot(table(sapply(mySeqs$seq, nchar))) ##Check seq length
  numChar<-sapply(mySeqs$seq, nchar) ##Record the length of each seq and store as numChar
  pullMe1<-grep(Y, numChar) ##Pull the indices of every seq length that is equal to 195
  pullMeAll<-c(pullMe1,pullMe2,pullMe3) ##Combine the two groups of indices
  pullNames<-names(fastafile)[pullMeAll] ##Pull the names of the previously pulled indices and place them in pullNames
  pullNames
  QC1file<-fastafile[c(which(names(fastafile) %in% pullNames))] ##Pulls seqs from original file for the seqs which names are found in the pullNames vector
  write.fasta(as.list(QC1file), names=myNames, file= paste0(Z, V, ".fasta")) ##Writes a new fasta file containing the seqs with correct length of NA
}
QC("Gill112017_HRSVA_G.pep", 297, "HRSVA_G", "/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/SBM/A/G/")
#############################################################################################################################################################################################################################################################
##Alignment
#############################################################################################################################################################################################################################################################
setwd("/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/SBM/A/G/")
Align<- function(X, Y){
  library(msa)
  orf <- readAAStringSet(file = X, format = "fasta") ##set name of fasta file containing the seqs with correct length/ previously translated to AA
  timestamp()
  musAlign<-msa(orf, "ClustalOmega") ##perform msa to align sequences
  timestamp()
  writeXStringSet(unmasked(musAlign), file= Y) ##set name of file
}
Align("HRSVA_G.fasta", "HRSVA_G_AA.phylip")
#############################################################################################################################################################################################################################################################
##Sequnce Based Mapping
#############################################################################################################################################################################################################################################################
SBM<-function(X, Y, Z){
  library(stringdist)
  viruAln<-read.alignment(file = X, format = "fasta", forceToLower = FALSE)
  table(sapply(viruAln$seq, nchar))
  plot(table(sapply(viruAln$seq, nchar)))
  for (i in 1:length(viruAln$nam)) {
    name<-viruAln$nam[[i]]
    seq<-viruAln$seq[[i]]
    assign(name, seq)
    #A epitope
    a<-substr(seq, 1, Y)
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
  write.csv(myMat, file= paste0(Z, "_DM.csv"))
  ##Dimension reduction to allow for 2D plotting
  ##Classical multidimensional scaling on distance matrix
  myMDS<-cmdscale(myMat, k=2, eig=T)
  ##Check goodness of fit
  myMDS$GOF
  ##Plot the reduced data
  pdf(file = paste0(Z, "SBM_plot.pdf"))
  plot(myMDS$points[,1],myMDS$points[,2])
  #text(myMDS$points[,1],myMDS$points[,2], labels=rownames(myMDS$points), cex = .5)
  dev.off()

}
setwd("/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/SBM/A/G")
SBM("HRSVA_G_AA.phylip", 322, "test_Gill_RSVA_G")
#############################################################################################################################################################################################################################################################
##Trait Association
#############################################################################################################################################################################################################################################################
myPheno<-read.csv("/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/myPheno.edited.csv", row.names = 1)
setwd("/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/SBM/DM_edited")
TA<- function(X, Y, Z){
  ##X is a data frame that contains the trait we are measuring association of
  ##Y is the name of the output stats .csv
  ##Z is the working directory of where PNG files should be written to
  library(Rtsne)
  library(vegan)
  myPheno$Global.Unique.Id<-as.character(myPheno$Global.Unique.Id) ##should be character vector, sets Name/ID col as character
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
    png(filename = paste0(Z, dir()[i], "_AnosimPlot.png"), width = 7, height = 7, units = "in", pointsize = 12, bg = "white", res = 1000, type = "cairo")
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
    png(filename = paste0(Z, dir()[i], "_Hallrsv_bDisper_sdEllipse.png"), width = 7, height = 7, units = "in", pointsize = 12, bg = "white", res = 1000, type = "cairo")
    try(plot(dispersion, hull=FALSE, ellipse=TRUE)) ##sd ellipse
    dev.off()
    myStats[i,1:4]<-c(myAnosim.pval, myAdonis.pval, bDis.perm.pval, bDis.anova.pval)
    myStats[i,5]<-dir()[i]
    print(dir()[i])
  }
  write.csv(myStats, paste0(Z, "TA_Stats_", Y, "_.csv"))
  return(myStats)
}
STATS<-TA(myPheno, "Aspires", "/Users/wbender/Library/CloudStorage/Box-Box/AndersonLabShared/William_Bender/RSV_Variation_Severity-Aspires/Aspires_Seqs/SBM/")

