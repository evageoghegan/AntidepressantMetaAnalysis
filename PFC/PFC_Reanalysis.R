library(devtools)
library(plyr)
library(tidyverse)
library(gemma.R)
library(dplyr)




setwd("/Users/evageoghegan/PFC")





setwd("/Users/evageoghegan/PFC")
library(dplyr)
DE_Results <- read.csv("GSE26836_PFC_Reanalysis.csv")
#The structure of the new DE_Results object:
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)

#Skip to step 9.9
GSE_ID<-"GSE26836"
ComparisonsOfInterest<-c("GSE26836_Amitriptyline")
NamesOfFoldChangeColumns<-c("Coef.treatment_factorAMI")
NamesOfTstatColumns<-c("t.treatment_factorAMI")
source("Function_CollapsingDEResults_OneResultPerGene.R")
#Example usage:

source("/Users/evageoghegan/PFC/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)


#Repeat for GSE84183
setwd("/Users/evageoghegan/PFC")
library(dplyr)
DE_Results <- read.csv("GSE84183_PFC_Reanalysis.csv")
#The structure of the new DE_Results object:
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)

#Skip to step 9.9
GSE_ID<-"GSE84183"
ComparisonsOfInterest<-c("GSE84183_fluoxetine")
NamesOfFoldChangeColumns<-c("Coef.treatment_factorFLX")
NamesOfTstatColumns<-c("t.treatment_factorFLX")
source("Function_CollapsingDEResults_OneResultPerGene.R")
#Example usage:
source("/Users/evageoghegan/PFC/Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#Reanalysis finished
