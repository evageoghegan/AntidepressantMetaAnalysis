#Adpating this code for Phi's results to only include the traditional antidepressant datasets


library(gemma.R)
library(tidyr)
library(dplyr)

getwd()
setwd("/Users/evageoghegan/ComparisonforPhi")
list.files()


#By hand, sort through the ResultSets_toScreen.csv to identify the statistical contrasts useful to your meta-analysis

#Delete all of the rows that you determine to be not applicable to your meta-analysis

#Save the file as a comma separate variable file named "ResultSets_Screened.csv"

#Afterwards, it would be good to double check your selections with me before proceeding forward.

#Done that, and a few have been excluded. 

#Saved in files as ResultSets_Screened.csv or excel version also

#Now read in that file with the function 

#For this code, have taken out all contrast IDs without fluoxetine/other trad antidepressants

ResultSet_contrasts<-read.csv("Result_SetsScreenedFluoxetine.csv", header=TRUE, stringsAsFactors = FALSE )
str(ResultSet_contrasts)

source("Function_DownloadingDEResults.R")

#example usage:
DownloadingDEResults(ResultSet_contrasts)




str(differentials)
#I have a list of 8
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
ComparisonsOfInterest<-c("GSE205325_fluoxetine")

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
ComparisonsOfInterest<-c("GSE56028_imipramine", "GSE56028_fluoxetine")

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
ComparisonsOfInterest<-c("GSE27532_desipramine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#5
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
ComparisonsOfInterest<-c("GSE63469_venlafaxine_high", "GSE63469_venlafaxine_low")
#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)


#6
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
ComparisonsOfInterest<-c("GSE81672_imipramine")

#Collapse into one result per gene
source("Function_CollapsingDEResults_OneResultPerGene.R")
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

##############
#Now to realign these results 
#Input mouse and rat datasets

setwd("/Users/evageoghegan/ComparisonforPhi")

############
library(plyr)
library(dplyr)
#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other
#HAs to be both mouse and rat sets


ListOfMouseDEResults<-list(DEResults_GSE118670, DEResults_GSE26836, DEResults_GSE84183,DEResults_GSE27532, DEResults_GSE63469, DEResults_GSE81672, DEResults_GSE43261_Ventral, DEResults_GSE43261_Dorsal)

AligningMouseDatasets(ListOfMouseDEResults)

#################

ListOfRatDEResults<-list(DEResults_GSE109445, DEResults_GSE205325, DEResults_GSE56028)

AligningRatDatasets(ListOfRatDEResults)


############

#Code for aligning the rat and mice results:

#This code isn't nicely functionalized yet
#It also assumes that there are mouse datasets
#It will break if there are only rat datasets - this needs to be fixed.

################

#First: What are gene orthologs?

# Homology refers to biological features including genes and their products that are descended from a feature present in a common ancestor.

# Homologous genes become separated in evolution in two different ways: separation of two populations with the ancestral gene into two species or gene duplication of the ancestral gene within a lineage:

### Genes separated by speciation are called orthologs.
### Genes separated by gene duplication events are called paralogs.

#This definition came from NCBI (https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html)


#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"

MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	28101 obs. of  7 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	28101 obs. of  7 variables:



#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")

MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")


#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:
colnames(MetaAnalysis_FoldChanges)

# [1] "Rat_EntrezGene.ID"                 "Mouse_EntrezGene.ID"              
# [3] "MouseVsRat_EntrezGene.ID"          "GSE126678_LPS_Acute"              
# [5] "GSE126678_LPS_SubchronicPlusAcute" "GSE126678_LPS_Subchronic"         
# [7] "GSE181285_LPS_Acute"               "GSE205325_LPS_Chronic" 

#Example scatter plot comparing two datasets:
plot(MetaAnalysis_FoldChanges$GSE81672_imipramine~MetaAnalysis_FoldChanges$GSE205325_fluoxetine)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))


## Meta Analysis
#This code caculates the number of NAs (i.e., the number of statistical contrasts lacking real differential expression results) in each row (i.e., for each gene):
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

#Or, since there are a limited number of integer answers, I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)

#For this dataset, I'm going to try running a meta-analysis using genes that were found in at least 4 sets of differential expression results
#Since there are 5 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#2 NA is too many
#NumberOfComparisons=5
#CutOffForNAs=2
#I have 5 statistical contrasts total (comparisons)
#and 2 NA is too many
NumberOfComparisons=13
CutOffForNAs=2
#set higher to 9 or 10
#cutoff for NAs lower for just SSRIs meta analysis -> 1 or 2


########################

#Running a basic meta-analysis:

#This function  is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV

source("Function_RunBasicMetaAnalysis.R")                         
#Example usage:

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)

colnames(MetaAnalysis_FoldChanges)
#once you find out which columns need to be removed, you can do it via number
#new object e.g metaanalysis fold changes SSRIs
#MetaAnalysis_FoldChanges_SSRIs<-MetaAnalysis_FoldChanges[,-c(3, 5, 8, 10)]
#the negative means "not those"
#do that for sampling variances and fold changes 
#MetaAnalysis_SV_SSRIs<-MetaAnalysis_SV[,-c(3, 5, 8, 10)]



#this will give you column names and what order theya re in
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...


#Example output:

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5 
# 13355  3200  5059   277  5293   917 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	16555 obs. of  8 variables:
#   $ Rat_EntrezGene.ID                 : chr  "114087" "191569" "246307" "65041" ...
# $ Mouse_EntrezGene.ID               : chr  "23825" "18585" "66514" "20480" ...
# $ MouseVsRat_EntrezGene.ID          : chr  "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# $ LPS_SubchronicAndAcute_vs_Vehicle : num  -0.0239 -0.0858 -0.0686 0.0891 0.0376 ...
# $ LPS_Acute_vs_Vehicle              : num  -0.0938 0.2524 -0.109 0.0788 -0.0475 ...
# $ LPS_Subchronic_vs_Vehicle         : num  0.0355 0.2735 0.0415 -0.0116 0.0299 ...
# $ LPS_Acute850ugPerKg_vs_Vehicle    : num  0.1754 0.1651 0.1878 -0.0146 0.103 ...
# $ LPS_Chronic_vs_Vehicle_AllStressed: num  0.223 -0.3824 -0.197 -0.0051 0.0393 ...
# NULL

str(metaOutput)
# num [1:16555, 1:6] 0.0234 0.1315 -0.023 0.0437 0.0769 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16555] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate         SE       pval        CI_lb      CI_ub Number_Of_Comparisons
# 23825_114087      0.02335239 0.04876404 0.63202016 -0.072223377 0.11892815                     5
# 18585_191569      0.13152931 0.06603932 0.04640599  0.002094618 0.26096399                     5
# 66514_246307     -0.02304066 0.05315084 0.66465470 -0.127214402 0.08113308                     5
# 20480_65041       0.04374446 0.03359984 0.19294219 -0.022110021 0.10959893                     5
# 13726_25437       0.07691445 0.04407249 0.08095348 -0.009466050 0.16329495                     5
# 16952_25380       0.04465958 0.10253423 0.66315765 -0.156303828 0.24562298                     5                   

tail(metaOutput)
# Log2FC_estimate         SE       pval       CI_lb      CI_ub Number_Of_Comparisons
# 108168734_NA      0.11566278 0.07742137 0.13519162 -0.03608031 0.26740587                     4
# 108168883_NA      0.16641238 0.09078313 0.06679127 -0.01151929 0.34434405                     4
# 108168923_NA      0.07757892 0.16987291 0.64789533 -0.25536587 0.41052370                     4
# 108168987_NA     -0.13375235 0.07407534 0.07097679 -0.27893734 0.01143264                     4
# 108169023_NA     -0.06106329 0.04794113 0.20276476 -0.15502618 0.03289960                     4
# 113002583_NA      0.02982367 0.07843726 0.70377969 -0.12391053 0.18355788                     4

########################################

#Accessing some additional gene annotation:

#Reading in a database containing more detailed gene annotation:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat_20240425.csv", header = TRUE, row.names = 1)

colnames(HOM_MouseVsRat)

#Renaming the columns so that we can easily join the annotation to our meta-analysis results:
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)

HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

#################

## Multiple comparison corrections

#This code runs a function that corrects the meta-analysis output to take into account the fact that we are running the statistical calculations thousands of times and therefore have a heightened risk of false discovery (false discovery rate correction) 

source("Function_FalseDiscoveryCorrection.R")

#Example usage:

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)



############################

#Next we will determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#Here are the top results as listed by mouse gene symbol:
metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:20)]
# "Fibcd1"   "Parm1"    "Cables1"  "Pla2g5"   "Gabrr2"   "Unc13c"   "Akap8l"   "Rgs14"   
# "Grm8"     "Slc4a5"   "Spmip7"   "Iqcg"     "Zfp691"   "Themis2"  "Aard"     "Tmem150b"
# "Tpbg"     "Chst10"   "Dtwd1"    "Iqch"    

#Or as listed by mouse entrez id:
metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:20)]
# "98970"  "231440" "63955"  "18784"  "14409"  "208898" "54194"  "51791"  "14823" 
# "232156" "73862"  "69707"  "195522" "230787" "239435" "330460" "21983"  "98388" 
# "69185"  "78250" 


#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:


#Range is mostly -1 to 1, but there are a few with Log2FC as big as -2-3

source("Function_MakeForestPlots.R")
#Note - this function currently uses Entrez ID (NCBI ID) as it's input
#It needs to be formatted as a character (not an integer) to work
#It can accept either mouse annotation (Entrez ID) or rat annotation (Entrez ID) as its input

#Example Usage:

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="208898", species="Mouse") 



MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="17289", species="Mouse") 










