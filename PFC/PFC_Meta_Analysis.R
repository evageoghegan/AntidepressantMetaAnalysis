#Eva rerun of Sophia's meta-analysis




#GSE84183 - reanalyzed, in your hippocampal folder
#GSE26836 - reanalyzed, in your hippocampal folder
#GSE230149 - pull from Gemma (just excluded before because it is parietal)

library(devtools)
library(plyr)
library(tidyverse)
library(gemma.R)
library(dplyr)




setwd("/Users/evageoghegan/PFC")



##############
#Running the Meta-analysis

#This code caculates the number of NAs (i.e., the number of statistical contrasts lacking real differential expression results) in each row (i.e., for each gene):
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

table(MetaAnalysis_FoldChanges_NAsPerRow)
#I have 17 contrasts
NumberOfComparisons=17
CutOffForNAs=5

##########
#Running a basic meta-analysis:
source("Function_RunBasicMetaAnalysis.R")      

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  

str(metaOutput)
head(metaOutput)
tail(metaOutput)

###############
#Accessing some additional gene annotation:

#Reading in a database containing more detailed gene annotation:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat_20240425.csv", header = TRUE, row.names = 1)

colnames(HOM_MouseVsRat)
#Renaming the columns so that we can easily join the annotation to our meta-analysis results:
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)
HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

##################
# Multiple comparison corrections
#This code runs a function that corrects the meta-analysis output to take into account the fact that we are running the statistical calculations thousands of times and therefore have a heightened risk of false discovery (false discovery rate correction) 

source("Function_FalseDiscoveryCorrection.R")

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)





#Make sure that multtest and plyr are loaded or you may come across an error
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")
library(plyr)

#Next we will determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#Here are the top results as listed by mouse gene symbol:
metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:20)]

#Or as listed by mouse entrez id:
metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:20)]

#Let's plot some of those top results!
hist(metaOutputFDR[,1], breaks=40)

setwd("/Users/evageoghegan/PFC")
library(metafor)
source("Function_MakeForestPlots.R")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="11966", species="Mouse") 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="20190", species="Mouse") 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="26557", species="Mouse") 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="16768", species="Mouse") 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="14824", species="Mouse") 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="232679", species="Mouse")

save.image("PFC_Meta_analysis.RData")












