library(gemma.R)
library(tidyr)
library(dplyr)

getwd()
setwd("/Users/evageoghegan/Brain Data Alchemy Project")
list.files()




MyDatasets_Screened<-read.csv("DatasetScreeningFinal.csv", stringsAsFactors = FALSE, header=TRUE)

str(MyDatasets_Screened)

#This allows us to see our screened datasets

#Ask if this is meant to be a csv file of the statistical contrast info or the original information about the dataset
#DatasetScreeningFinal is the list of 17 experiments that I ran statistical contrasts for 

#Now pull out experiment IDs:
ExperimentIDs<-MyDatasets_Screened$experiment.shortName
ExperimentIDs

#This code makes a data.frame that includes all of the contrast ids for each dataset with their basic metadata in a format that is easily readable in a spreadsheet program
#This is how when you look at all the information from the datasets to review what it is being subsetted by

source("Function_GettingResultSetInfoForDatasets.R")

GettingResultSetInfoForDatasets(ExperimentIDs)
#By hand, sort through the ResultSets_toScreen.csv to identify the statistical contrasts useful to your meta-analysis

#Delete all of the rows that you determine to be not applicable to your meta-analysis

#Save the file as a comma separate variable file named "ResultSets_Screened.csv"

#Afterwards, it would be good to double check your selections with me before proceeding forward.

#Done that, and a few have been excluded. 

#Saved in files as ResultSets_Screened.csv or excel version also

#Now read in that file with the function 

ResultSet_contrasts<-read.csv("ResultSets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )
str(ResultSet_contrasts)

source("Function_DownloadingDEResults.R")

#example usage:
DownloadingDEResults(ResultSet_contrasts)


#This has been updated to include tianeptine (originally missed because named differently in stat constrats) 18th #November 2023

str(differentials)
#I have a list of 17
#Each of those is a result set containing all of the differential expression results for that particular variable. 
#If the variable has more than one level (e.g., 3 LPS dosages), you will have several statistical contrasts included in the result set differential expression output.  
#Here is how you can access and review the differential expression results for a particular result set id:
str(differentials[1])
#For record-keeping purposes, let's save the differential expression results for each result set:

#The reason for this is because we will have a hard time replicating our results in the future with just our code
#Because Gemma updates the annotation on the data every time a new reference genome is released

source("Function_SavingGemmaDEResults_forEachResultSet.R")

#Example usage:

SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
#Next we will start working with cleaning up the results for a single result set

#We need to filter down our differential expression results to just the rows with good gene annotation

DE_Results<-differentials[[1]]

str(DE_Results)

source("Function_FilteringDEResults_GoodAnnotation.R")

#Example of using the function for a dataset:

FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)

#Next we are going to pull out the differential expression for the specific statistical contrasts that we are interested in

source("Function_ExtractingDEResultsForContrasts.R")

#Example usage:
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE109445_fluoxetine")

#####################

#Sometimes gene expression is measured using multiple probes (microarray)
#Next we need to collapse our differential expression results down to one result per gene
#At the same time, we will calculate the standard error for our effect size (Log2FC) using the t-statistic
#And then square the standard error to get the sampling variance

source("Function_CollapsingDEResults_OneResultPerGene.R")

#Notes about parameters for function CollapsingDEResults_OneResultPerGene()
#GSE_ID is a string indicating the name of the Gemma dataset
#DE_Results_GoodAnnotation is the data frame outputted by our previous function
#ComparisonsOfInterest is a character vector containing the names of the group comparisons of interest within this dataset. Important: These group comparisons should be listed in exactly the same order as the order that you provide the column names for their associated Fold Change and Tstat output.
#NamesOfFoldChangeColumns is a vector containing the names of the columns of DE_Results_GoodAnnotation containing the FoldChange results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector
#NamesOfTstatColumns is a vector containing the names of the columns of DE_Results_GoodAnnotation containing the Tstat results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector

#Example usage:
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

# [1] "Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 21614
# [1] "# of rows with unique Gene Symbols:"
# [1] 21614
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 21614     3
# [1] "Output: Named DEResults_GSE126678"

#################


###Then you can move onto the other differentials.

##Number 2
str(differentials[2])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[2]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE123027_ECT")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

##Number 3
str(differentials[3])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[3]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE205325_fluoxetine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

##Number 4
str(differentials[4])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[4]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE230148_TMS_modifier_theta_burst", "GSE230148_TMS_1Hz")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

##Number 5
str(differentials[5])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[5]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE_230149_TMS_intermittent_theta_burst")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[6])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[6]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE27532_desipramine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[7])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[7]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE56028_imipramine","GSE56028_tianeptine", "GSE56028_fluoxetine", "GSE56028_agomelatine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[8])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[8]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE61301_imipramine_yohimbine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[9])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[9]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE63469_venlafaxine_high", "GSE63469_venlafaxine_low")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[10])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[10]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE73798_ketamine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

str(differentials[11])
source("Function_SavingGemmaDEResults_forEachResultSet.R")
SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
DE_Results<-differentials[[11]]
str(DE_Results)
source("Function_FilteringDEResults_GoodAnnotation.R")
FilteringDEResults_GoodAnnotation(DE_Results)
str(DE_Results_GoodAnnotation)
source("Function_ExtractingDEResultsForContrasts.R")
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Rename the results of interest 
ComparisonsOfInterest<-c("GSE81672_ketamine", "GSE81672_imipramine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)










