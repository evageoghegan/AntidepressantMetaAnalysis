#1) combining the gene annotation file and limma results for each hippocampal dataset into a single file
#2) save it as a .csv
#3) read it into R using the read.csv() function

setwd("/Users/evageoghegan/Brain Data Alchemy Project")
library(dplyr)
DE_Results <- read.csv("GSE118670.csv")
#The structure of the new DE_Results object:
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)

#Skip to step 9.9
GSE_ID<-"GSE118670"
ComparisonsOfInterest<-c("GSE118670_Fluoxetine")
NamesOfFoldChangeColumns<-c("Coef.treatment_factorFLX")
NamesOfTstatColumns<-c("t.treatment_factorFLX")
source("Function_CollapsingDEResults_OneResultPerGene.R")
#Example usage:

source("/Users/evageoghegan/Brain Data Alchemy Project/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)


#Repeat for GSE26836
setwd("/Users/evageoghegan/Brain Data Alchemy Project")
library(dplyr)
DE_Results <- read.csv("GSE26836.csv")
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
source("/Users/evageoghegan/Brain Data Alchemy Project/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)


#Repeat for GSE84183
setwd("/Users/evageoghegan/Brain Data Alchemy Project")
library(dplyr)
DE_Results <- read.csv("GSE84183.csv")
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
source("/Users/evageoghegan/Brain Data Alchemy Project/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#Repeat for GSE43261_Ventral
setwd("/Users/evageoghegan/Brain Data Alchemy Project")
library(dplyr)
DE_Results <- read.csv("GSE43261_Ventral.csv")
#The structure of the new DE_Results object:
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
#had done up until here
#Skip to step 9.9
GSE_ID<-"GSE43261_Ventral"
ComparisonsOfInterest<-c("GSE43261_Ventral_fluoxetine")
NamesOfFoldChangeColumns<-c("Coef.treatment_factorSSRI")
NamesOfTstatColumns<-c("t.treatment_factorSSRI")
source("Function_CollapsingDEResults_OneResultPerGene.R")
#Example usage:
source("/Users/evageoghegan/Brain Data Alchemy Project/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#Repeat for GSE43261_Dorsal
setwd("/Users/evageoghegan/Brain Data Alchemy Project")
library(dplyr)
DE_Results <- read.csv("GSE43261_Dorsal.csv")
#The structure of the new DE_Results object:
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)

#Skip to step 9.9
GSE_ID<-"GSE43261_Dorsal"
ComparisonsOfInterest<-c("GSE43261_Dorsal_fluoxetine")
NamesOfFoldChangeColumns<-c("Coef.treatment_factorSSRI")
NamesOfTstatColumns<-c("t.treatment_factorSSRI")
source("Function_CollapsingDEResults_OneResultPerGene.R")
#Example usage:
source("/Users/evageoghegan/Brain Data Alchemy Project/Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)




citation("plyr")
citation("dplyr")
citation("gemma.R")
citation("limma")
citation("fgsea")
packageVersion("plyr")
packageVersion("dplyr")
packageVersion("gemma.R")
