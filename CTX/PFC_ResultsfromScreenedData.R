setwd("/Users/evageoghegan/PFC")

MyDatasets_Screened<-read.csv("PFC_MyDatasets_Screened.csv", stringsAsFactors = FALSE, header=TRUE)

str(MyDatasets_Screened)

#Read in your screened result sets:

ResultSet_contrasts<-read.csv("PFC_MyDatasets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )
str(ResultSet_contrasts)

source("Function_DownloadingDEResults.R")
DownloadingDEResults(ResultSet_contrasts)

str(differentials)
#Here is how you can access and review the differential expression results for a particular result set id:
str(differentials[1])

#For record-keeping purposes, let's save the differential expression results for each result set:

source("Function_SavingGemmaDEResults_forEachResultSet.R")

#######################

#Next we will start working with cleaning up the results for a single result set

DE_Results<-differentials[[1]]

str(DE_Results)

#Reading in the function

source("Function_FilteringDEResults_GoodAnnotation.R")

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

###################
#Next we are going to pull out the differential expression for the specific statistical contrasts that we are interested in

source("Function_ExtractingDEResultsForContrasts.R")

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE28644_fluoxetine")


####################
#Next we need to collapse our differential expression results down to one result per gene
#At the same time, we will calculate the standard error for our effect size (Log2FC) using the t-statistic
#And then square the standard error to get the sampling variance

source("Function_CollapsingDEResults_OneResultPerGene.R")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#Do this with the remaining contrasts:
DE_Results<-differentials[[2]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE93041_ketamine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[3]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE81672_ketamine", "GSE81672_imipramine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[4]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE150264_imipramine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[5]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE84185_fluoxetine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[6]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE168172_duloxetine", "GSE168172_sertraline")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[7]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE138802_ketamine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[8]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE129359_duloxetine")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[9]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE45229_quetiapine_low_dose", "GSE45229_quetiapine_high_dose")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[10]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE253280_MDMA")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

DE_Results<-differentials[[11]]
str(DE_Results)
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
ComparisonsOfInterest<-c("GSE230149_TMS")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)



