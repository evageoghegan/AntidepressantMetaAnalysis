#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include co-variates
#Updated version: March 10, 2026

######################


#Trying out a mondo-meta-analysis with both hpc and cortex 


######################

#Grabbing input data and setting the working directory:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/R_Code_And_Workspaces/Meta-analysis/HPC_November_Final_Meta_Analysis.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MondoMeta_PFCvsHPC")

######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

#For documentation:
write.csv(MetaAnalysis_FoldChanges, "MetaAnalysis_FoldChanges_HPC.csv")
write.csv(MetaAnalysis_FoldChanges_ForMeta, "MetaAnalysis_FoldChanges_ForMeta_HPC.csv")
write.csv(MetaAnalysis_SV, "MetaAnalysis_SV_HPC.csv")
write.csv(MetaAnalysis_SV_ForMeta, "MetaAnalysis_SV_ForMeta_HPC.csv")

MetaAnalysis_FoldChanges_HPC<-MetaAnalysis_FoldChanges
MetaAnalysis_SV_HPC<-MetaAnalysis_SV

colnames(MetaAnalysis_FoldChanges_HPC)[c(4:25)]<-paste(colnames(MetaAnalysis_FoldChanges_HPC)[c(4:25)], "HPC", sep="_")

colnames(MetaAnalysis_SV_HPC)[c(4:25)]<-paste(colnames(MetaAnalysis_SV_HPC)[c(4:25)], "HPC", sep="_")

HPC_Covariates<-read.csv("HPC_Covariates.csv", header=TRUE, stringsAsFactors = FALSE)

str(HPC_Covariates)
# 'data.frame':	22 obs. of  6 variables:
# $ X              : chr  "GSE123027_ECT" "GSE27532_desipramine" "GSE63469_venlafaxine_high" "GSE63469_venlafaxine_low" ...
# $ ADType         : num  0.5 -0.5 -0.5 -0.5 0.5 0.5 -0.5 -0.5 -0.5 -0.5 ...
# $ Dissection     : num  -0.5 0.5 -0.5 -0.5 -0.5 -0.5 -0.5 0.5 -0.5 0.5 ...
# $ Platform       : num  -0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ DepressionModel: num  0.5 0.5 0.5 0.5 -0.5 0.5 0.5 0.5 -0.5 -0.5 ...
# $ SampleSize     : int  24 16 4 4 24 17 18 32 6 16 ...

HPC_Covariates$Study_Treatment<-paste(HPC_Covariates$X, "HPC", sep="_")
HPC_Covariates$Study_Treatment

MetaAnalysis_FoldChanges_PFC<-read.csv("MetaAnalysis_FoldChanges_PFC.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(MetaAnalysis_FoldChanges_PFC)
MetaAnalysis_FoldChanges_PFC<-MetaAnalysis_FoldChanges_PFC[,-1]
colnames(MetaAnalysis_FoldChanges_PFC)[c(4:19)]<-paste(colnames(MetaAnalysis_FoldChanges_PFC)[c(4:19)], "PFC", sep="_")

MetaAnalysis_SV_PFC<-read.csv("MetaAnalysis_SV_PFC.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(MetaAnalysis_SV_PFC)
MetaAnalysis_SV_PFC<-MetaAnalysis_SV_PFC[,-1]
colnames(MetaAnalysis_SV_PFC)[c(4:19)]<-paste(colnames(MetaAnalysis_SV_PFC)[c(4:19)], "PFC", sep="_")

PFC_Covariates<-read.csv("PFC_Covariates.csv", header=TRUE, stringsAsFactors = FALSE)

str(PFC_Covariates)

# 'data.frame':	16 obs. of  7 variables:
#   $ Study_Treatment: chr  "GSE26836_Amitriptyline" "GSE84183_fluoxetine" "GSE118670_Fluoxetine" "GSE28644_fluoxetine" ...
# $ AD.Type        : num  -0.5 -0.5 -0.5 -0.5 0.5 0.5 -0.5 -0.5 -0.5 -0.5 ...
# $ Dissection_Acg : num  0 0.5 -0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 -0.5 ...
# $ Dissection_PFC : num  0 -0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5 0.5 0.5 ...
# $ Platform       : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 -0.5 -0.5 ...
# $ DepressionModel: num  -0.5 0.5 -0.5 0.5 -0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ SampleSize     : int  6 32 16 33 6 18 19 36 9 9 ...

PFC_Covariates$Study_Treatment<-paste(PFC_Covariates$Study_Treatment, "PFC", sep="_")
PFC_Covariates$Study_Treatment

Covariates<-data.frame(Study_Treatment=c(HPC_Covariates$Study_Treatment, PFC_Covariates$Study_Treatment), ADType=c(HPC_Covariates$ADType, PFC_Covariates$AD.Type), Dissection=c(rep(0.5, length(HPC_Covariates$Study_Treatment)), rep(-0.5, length(PFC_Covariates$Study_Treatment))),  Platform=c(HPC_Covariates$Platform, PFC_Covariates$Platform), DepressionModel=c(HPC_Covariates$DepressionModel, PFC_Covariates$DepressionModel))

colnames(MetaAnalysis_FoldChanges_PFC)

MetaAnalysis_FoldChanges<-join(MetaAnalysis_FoldChanges_HPC, MetaAnalysis_FoldChanges_PFC, by="MouseVsRat_EntrezGene.ID", type="inner")

str(MetaAnalysis_FoldChanges)
colnames(MetaAnalysis_FoldChanges) #Rat and mouse entrez id are in there twice
MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[,-c(26:27)]
colnames(MetaAnalysis_FoldChanges)
# [1] "Rat_EntrezGene.ID"                           "Mouse_EntrezGene.ID"                        
# [3] "MouseVsRat_EntrezGene.ID"                    "GSE123027_ECT_HPC"                          
# [5] "GSE27532_desipramine_HPC"                    "GSE63469_venlafaxine_high_HPC"              
# [7] "GSE63469_venlafaxine_low_HPC"                "GSE73798_ketamine_HPC"                      
# [9] "GSE81672_ketamine_HPC"                       "GSE81672_imipramine_HPC"                    
# [11] "GSE84183_fluoxetine_HPC"                     "GSE26836_Amitriptyline_HPC"                 
# [13] "GSE118670_Fluoxetine_HPC"                    "GSE43261_Dorsal_fluoxetine_HPC"             
# [15] "GSE43261_Ventral_fluoxetine_HPC"             "GSE109445_fluoxetine_HPC"                   
# [17] "GSE205325_fluoxetine_HPC"                    "GSE_230149_TMS_intermittent_theta_burst_HPC"
# [19] "GSE230148_TMS_modifier_theta_burst_HPC"      "GSE230148_TMS_1Hz_HPC"                      
# [21] "GSE56028_imipramine_HPC"                     "GSE56028_tianeptine_HPC"                    
# [23] "GSE56028_fluoxetine_HPC"                     "GSE56028_agomelatine_HPC"                   
# [25] "GSE61301_imipramine_yohimbine_HPC"           "GSE26836_Amitriptyline_PFC"                 
# [27] "GSE84183_fluoxetine_PFC"                     "GSE118670_Fluoxetine_PFC"                   
# [29] "GSE28644_fluoxetine_PFC"                     "GSE93041_ketamine_PFC"                      
# [31] "GSE81672_ketamine_PFC"                       "GSE81672_imipramine_PFC"                    
# [33] "GSE150264_imipramine_PFC"                    "GSE168172_duloxetine_PFC"                   
# [35] "GSE168172_sertraline_PFC"                    "GSE138802_ketamine_PFC"                     
# [37] "GSE129359_duloxetine_PFC"                    "GSE45229_quetiapine_low_dose_PFC"           
# [39] "GSE45229_quetiapine_high_dose_PFC"           "GSE230149_TMS_PFC"                          
# [41] "GSE253280_MDMA_PFC"      

str(MetaAnalysis_FoldChanges)
#'data.frame':	30240 obs. of  41 variables:

MetaAnalysis_SV<-join(MetaAnalysis_SV_HPC, MetaAnalysis_SV_PFC, by="MouseVsRat_EntrezGene.ID", type="inner")

str(MetaAnalysis_SV)
colnames(MetaAnalysis_SV) #Rat and mouse entrez id are in there twice
MetaAnalysis_SV<-MetaAnalysis_SV[,-c(26:27)]
colnames(MetaAnalysis_SV)

# [1] "Rat_EntrezGene.ID"                           "Mouse_EntrezGene.ID"                        
# [3] "MouseVsRat_EntrezGene.ID"                    "GSE123027_ECT_HPC"                          
# [5] "GSE27532_desipramine_HPC"                    "GSE63469_venlafaxine_high_HPC"              
# [7] "GSE63469_venlafaxine_low_HPC"                "GSE73798_ketamine_HPC"                      
# [9] "GSE81672_ketamine_HPC"                       "GSE81672_imipramine_HPC"                    
# [11] "GSE84183_fluoxetine_HPC"                     "GSE26836_Amitriptyline_HPC"                 
# [13] "GSE118670_Fluoxetine_HPC"                    "GSE43261_Dorsal_fluoxetine_HPC"             
# [15] "GSE43261_Ventral_fluoxetine_HPC"             "GSE109445_fluoxetine_HPC"                   
# [17] "GSE205325_fluoxetine_HPC"                    "GSE_230149_TMS_intermittent_theta_burst_HPC"
# [19] "GSE230148_TMS_modifier_theta_burst_HPC"      "GSE230148_TMS_1Hz_HPC"                      
# [21] "GSE56028_imipramine_HPC"                     "GSE56028_tianeptine_HPC"                    
# [23] "GSE56028_fluoxetine_HPC"                     "GSE56028_agomelatine_HPC"                   
# [25] "GSE61301_imipramine_yohimbine_HPC"           "GSE26836_Amitriptyline_PFC"                 
# [27] "GSE84183_fluoxetine_PFC"                     "GSE118670_Fluoxetine_PFC"                   
# [29] "GSE28644_fluoxetine_PFC"                     "GSE93041_ketamine_PFC"                      
# [31] "GSE81672_ketamine_PFC"                       "GSE81672_imipramine_PFC"                    
# [33] "GSE150264_imipramine_PFC"                    "GSE168172_duloxetine_PFC"                   
# [35] "GSE168172_sertraline_PFC"                    "GSE138802_ketamine_PFC"                     
# [37] "GSE129359_duloxetine_PFC"                    "GSE45229_quetiapine_low_dose_PFC"           
# [39] "GSE45229_quetiapine_high_dose_PFC"           "GSE230149_TMS_PFC"                          
# [41] "GSE253280_MDMA_PFC" 

str(MetaAnalysis_SV)
#'data.frame':	30240 obs. of  41 variables

str(Covariates)
# 'data.frame':	38 obs. of  5 variables:
#   $ Study_Treatment: chr  "GSE123027_ECT_HPC" "GSE27532_desipramine_HPC" "GSE63469_venlafaxine_high_HPC" "GSE63469_venlafaxine_low_HPC" ...
# $ ADType         : num  0.5 -0.5 -0.5 -0.5 0.5 0.5 -0.5 -0.5 -0.5 -0.5 ...
# $ Dissection     : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ Platform       : num  -0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ DepressionModel: num  0.5 0.5 0.5 0.5 -0.5 0.5 0.5 0.5 -0.5 -0.5 ...

write.csv(Covariates, "Covariates_HPC_PFC.csv")
write.csv(MetaAnalysis_FoldChanges, "MetaAnalysis_FoldChanges_HPC_PFC_Shared.csv")
write.csv(MetaAnalysis_SV, "MetaAnalysis_SV_HPC_PFC_Shared.csv")

ADType<-Covariates$ADType
Dissection<-Covariates$Dissection
Platform<-Covariates$Platform
DepressionModel<-Covariates$DepressionModel

#Platform may make this crash - many of these genes aren't found in the old microarray datasets

Dissection
# [1]  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5
# [19]  0.5  0.5  0.5  0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5
# [37] -0.5 -0.5

######################

#Grabbing input data and setting the working directory:

# load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MondoMetaPFCvsHPC_DetailedAD/Workspace_MetaRegression_HPC_CTX_ADDetailed.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_HPCvsCTX")


#########

#Function:

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  #The function first provides information about how many of the statistical contrasts have NA values as their differential expression results for each gene:
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  #Then any row (gene) that has too many NAs is removed from the analysis:
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  #This matrix was originally 6 columns
  #It was made larger to incorporate heterogeneity statistics (8 stats)
  #And then even larger to include stats for the covariates
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 22)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    print(i)
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~Dissection, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~Dissection, vi=var)
      metaOutput[i, c(1:2)]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, c(3:4)]<-TempMeta$se #gives standard error
      metaOutput[i, c(5:6)]<-TempMeta$pval #gives pval
      metaOutput[i, c(7:8)]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, c(9:10)]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 11]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      metaOutput[i, 12]<-TempMeta$k #Metafor output: number of studies (contrasts) included in the analysis - sanity check, should be the same as column 6
      metaOutput[i, 13]<-TempMeta$p #Metafor output: number of coefficients in model
      metaOutput[i, 14]<-TempMeta$tau2 #estimated amount of (residual) heterogeneity
      metaOutput[i, 15]<-TempMeta$se.tau2 #SE of the estimated amount of (residual) heterogeneity
      metaOutput[i, 16]<-TempMeta$QE #When moderators are included in the model, this is the QE-test for residual heterogeneity
      metaOutput[i, 17]<-TempMeta$QEp #p-value of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 18]<-TempMeta$I2 #the I 2 statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes can be attributed to heterogeneity among the true effects
      metaOutput[i, 19]<-TempMeta$H2 #the H2 statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability
      metaOutput[i, 20]<-TempMeta$QM #test statistic of the omnibus test of moderators
      metaOutput[i, 21]<-TempMeta$QMp #corresponding p-value of the omnibus test of moderators
      metaOutput[i, 22]<-TempMeta$R2 #amount of heterogeneity accounted for by the moderators
      
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_AD_vs_Ctrl", "Log2FC_Dissection_HPCvsPFC",          
                          "SE_AD", "SE_Dissection",  
                          "pval_AD", "pval_Dissection",
                          "CI_lb_AD", "CI_lb_Dissection",
                          "CI_ub_AD", "CI_ub_Dissection",
                          "Number_Of_Comparisons", "Number_of_ Contrasts", "Number_of_Coefficients", "tau2_ResidualHeterogeneity", "SE_tau2_ResidualHeterogeneity", "QE_CochransQ_Teststat", "QEp_CochransQ_pval", "I2_PercentVar_TrueHeterogeneity", "H2_Ratio_EffectHetero_overSamplVar", "QM_ModeratorOmnibusTest", "QMp_Pval_ModeratorOmnibusTest", "R2_ModeratorVarianceExplained")
  
  #The row names for our output are the combined mouse-rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  #... and provide the user with an update about the newly created object:
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]

  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
}

######################

#Example Usage:
print(table(MetaAnalysis_FoldChanges_NAsPerRow))
MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 5972  620 2455 1083  811  372  316  205  134  889  226  554  790  312  472  645  398  342  351 
# 19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37 
# 223  163  182  178  207  678  537  996  695 1572 1411  644  578  414  771 1002 2270  732   22 
# 38 

NumberOfComparisons=38
CutOffForNAs=10
#I have 38 comparisons
#I want to make sure that I have studies from both HPC (22 comparisons) and CTX (16 comparisons)
#The cut off used for the CTX was 6 NA (so that there were at least 11 comparisons)
#That seems a little too strict for this analysis:
5972+620+2455+1083+811+372
#[1] 11313
#So that would only leave 11,000 genes
#Using a cut-off of 10 would get us:
5972+620+2455+1083+811+372+316+205+134+889
#[1] 12857
#...that's a little bit better?

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

######################

#Example Output:

# [1] "metaOutput:"
# num [1:12857, 1:22] -0.003115 -0.052317 -0.030794 -0.001621 -0.000629 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:12857] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:22] "Log2FC_AD_vs_Ctrl" "Log2FC_Dissection_HPCvsPFC" "SE_AD" "SE_Dissection" ...
# NULL
# [1] "Top of metaOutput:"
# Log2FC_AD_vs_Ctrl Log2FC_Dissection_HPCvsPFC       SE_AD SE_Dissection   pval_AD
# 23825_114087      -0.003115431                0.008478798 0.012360491    0.02472098 0.8010043
# 18585_191569      -0.052316932               -0.002286042 0.033803033    0.06760607 0.1216947
# 66514_246307      -0.030793775               -0.039018432 0.019266255    0.03853251 0.1099702
# 20480_65041       -0.001621439               -0.052079057 0.012225227    0.02445045 0.8944855
# 13726_25437       -0.000629093                0.003424980 0.009267352    0.01853470 0.9458790
# 16952_25380        0.041891287               -0.093306578 0.036841472    0.07368294 0.2555095
# pval_Dissection    CI_lb_AD CI_lb_Dissection    CI_ub_AD CI_ub_Dissection
# 23825_114087      0.73161364 -0.02734155      -0.03997344 0.021110685      0.056931031
# 18585_191569      0.97302535 -0.11856966      -0.13479150 0.013935795      0.130219412
# 66514_246307      0.31124615 -0.06855494      -0.11454076 0.006967391      0.036503900
# 20480_65041       0.03317299 -0.02558244      -0.10000106 0.022339564     -0.004157049
# 13726_25437       0.85339576 -0.01879277      -0.03290237 0.017534583      0.039752333
# 16952_25380       0.20539663 -0.03031667      -0.23772249 0.114099245      0.051109338
# Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 23825_114087                    38                   38                      2
# 18585_191569                    36                   36                      2
# 66514_246307                    38                   38                      2
# 20480_65041                     38                   38                      2
# 13726_25437                     38                   38                      2
# 16952_25380                     36                   36                      2
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 23825_114087               0.0012894072                  0.0011575719              48.53904
# 18585_191569               0.0304018167                  0.0096140177             187.54109
# 66514_246307               0.0087219949                  0.0031316330             121.96765
# 20480_65041                0.0017847076                  0.0011977888              54.88551
# 13726_25437                0.0001085089                  0.0005620116              45.82420
# 16952_25380                0.0166957432                  0.0084836881              66.83078
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity
# 23825_114087       7.912920e-02                        25.58809
# 18585_191569       3.877490e-23                        83.96617
# 66514_246307       2.824742e-11                        74.10034
# 20480_65041        2.274582e-02                        36.13572
# 13726_25437        1.263888e-01                         3.22258
# 16952_25380        6.524876e-04                        53.76416
# H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 23825_114087                           1.343871             0.117635154
# 18585_191569                           6.236813             0.001143397
# 66514_246307                           3.861054             1.025380422
# 20480_65041                            1.565820             4.536828868
# 13726_25437                            1.033299             0.034146398
# 16952_25380                            2.162824             1.603579846
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 23825_114087                    0.73161364                       0.00000
# 18585_191569                    0.97302535                       0.00000
# 66514_246307                    0.31124615                       0.00000
# 20480_65041                     0.03317299                      23.27257
# 13726_25437                     0.85339576                       0.00000
# 16952_25380                     0.20539663                       0.00000
# [1] "Bottom of metaOutput"
# Log2FC_AD_vs_Ctrl Log2FC_Dissection_HPCvsPFC       SE_AD SE_Dissection    pval_AD
# 240066_362845        -0.042965333                0.025245468 0.024633857    0.04926771 0.08113158
# 71640_100233177      -0.009613696               -0.019826794 0.019311361    0.03862272 0.61860676
# 232853_308320         0.033249577               -0.083072885 0.023564851    0.04712970 0.15824996
# 101197_312301        -0.003547277                0.004598907 0.009403619    0.01880724 0.70600670
# 22776_308322         -0.021476309               -0.046443265 0.020542375    0.04108475 0.29580867
# 243308_304336        -0.002815173               -0.032365071 0.014023173    0.02804635 0.84089288
# pval_Dissection    CI_lb_AD CI_lb_Dissection    CI_ub_AD CI_ub_Dissection
# 240066_362845        0.60836127 -0.09124681      -0.07131748 0.005316139      0.121808414
# 71640_100233177      0.60770977 -0.04746327      -0.09552594 0.028235877      0.055872350
# 232853_308320        0.07796055 -0.01293668      -0.17544540 0.079435836      0.009299633
# 101197_312301        0.80682147 -0.02197803      -0.03226260 0.014883479      0.041460418
# 22776_308322         0.25829681 -0.06173862      -0.12696789 0.018786006      0.034081365
# 243308_304336        0.24850618 -0.03030009      -0.08733490 0.024669740      0.022604755
# Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 240066_362845                      36                   36                      2
# 71640_100233177                    31                   31                      2
# 232853_308320                      38                   38                      2
# 101197_312301                      38                   38                      2
# 22776_308322                       34                   34                      2
# 243308_304336                      38                   38                      2
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 240066_362845                  0.012267080                  0.0047084878             111.15957
# 71640_100233177                0.005182678                  0.0027479632              61.81899
# 232853_308320                  0.011578090                  0.0045822798             100.03178
# 101197_312301                  0.000000000                  0.0006187254              29.88786
# 22776_308322                   0.003606176                  0.0027260244              42.97341
# 243308_304336                  0.001239626                  0.0013695455              37.18558
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity
# 240066_362845         4.012039e-10                        71.35941
# 71640_100233177       3.648678e-04                        55.72241
# 232853_308320         6.114145e-08                        66.53471
# 101197_312301         7.535919e-01                         0.00000
# 22776_308322          9.314919e-02                        44.97439
# 243308_304336         4.142308e-01                        20.05099
# H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 240066_362845                             3.491548              0.26256815
# 71640_100233177                           2.258479              0.26352344
# 232853_308320                             2.988170              3.10691381
# 101197_312301                             1.000000              0.05979423
# 22776_308322                              1.817336              1.27786271
# 243308_304336                             1.250797              1.33168218
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 240066_362845                      0.60836127                      0.000000
# 71640_100233177                    0.60770977                      0.000000
# 232853_308320                      0.07796055                     10.940614
# 101197_312301                      0.80682147                      0.000000
# 22776_308322                       0.25829681                     15.625116
# 243308_304336                      0.24850618                      5.255998


write.csv(metaOutput, "metaOutput_wCovariates_HPCvsCTX.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wCovariates_HPCvsCTX.csv")

colnames(metaOutput)

# [1] "Log2FC_AD_vs_Ctrl"                  "Log2FC_Dissection_HPCvsPFC"        
# [3] "SE_AD"                              "SE_Dissection"                     
# [5] "pval_AD"                            "pval_Dissection"                   
# [7] "CI_lb_AD"                           "CI_lb_Dissection"                  
# [9] "CI_ub_AD"                           "CI_ub_Dissection"                  
# [11] "Number_Of_Comparisons"              "Number_of_ Contrasts"              
# [13] "Number_of_Coefficients"             "tau2_ResidualHeterogeneity"        
# [15] "SE_tau2_ResidualHeterogeneity"      "QE_CochransQ_Teststat"             
# [17] "QEp_CochransQ_pval"                 "I2_PercentVar_TrueHeterogeneity"   
# [19] "H2_Ratio_EffectHetero_overSamplVar" "QM_ModeratorOmnibusTest"           
# [21] "QMp_Pval_ModeratorOmnibusTest"      "R2_ModeratorVarianceExplained"  

#############################

#I should probably run an FDR correction for all of these pvalues
#That function needs to be adapted too


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")                                 
library(multtest)


dim(metaOutput)
colnames(metaOutput)

####################

#Function:

FalseDiscoveryCorrection<-function(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation){
  
  #For the meta-analysis pval:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,5], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[23]<-"AD_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the ADType:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,6], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[24]<-"Dissection_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the QEp:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,17], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[25]<-"QEp_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the QMp:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,21], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[26]<-"QEp_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  
  #These results are returned to our global environment:
  metaOutputFDR<<-metaOutputFDR
  
  #We let the user know the basic structure of the meta-analysis output with FDR added to it (just to make sure everything still looks good)
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  #Then we make a dataframe that adds the annotation to that output:
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  #And then adds even more detailed gene annotation:
  
  #First the detailed annotation for the mouse genes:
  TempDF2<-join(TempDF, HOM_MouseVsRat[,c(4:5,9:11)], by="Mouse_EntrezGene.ID", type="left", match="first")
  
  #Next the annnotation for the rat genes:
  TempDF3<-join(TempDF2, HOM_MouseVsRat[,c(15:16,20:22)], by="Rat_EntrezGene.ID", type="left", match="first")
  
  #This is renamed and returned to our global environment:
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  #And written out into our working directory:
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #Then we make a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  #And give the user some information about their results:
  
  print("Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?")
  print(sum(metaOutputFDR_annotated[,25]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,25]),]))
  
  
}


FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)


# [1] "metaOutputFDR:"
# num [1:12857, 1:26] -0.003115 -0.052317 -0.030794 -0.001621 -0.000629 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:12857] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:26] "Log2FC_AD_vs_Ctrl" "Log2FC_Dissection_HPCvsPFC" "SE_AD" "SE_Dissection" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 29
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_AD_vs_Ctrl Log2FC_Dissection_HPCvsPFC
# 81897_338457             338457               81897        0.07491842                0.042918663
# 209239_307893            307893              209239       -0.06384987                0.009848819
# 75590_293847             293847               75590        0.06325243                0.028802933
# 217666_314196            314196              217666       -0.04027109               -0.015623669
# 21952_29388               29388               21952        0.12850947                0.032129480
# 69893_298377             298377               69893        0.06481482                0.019851656
# SE_AD SE_Dissection      pval_AD pval_Dissection    CI_lb_AD CI_lb_Dissection
# 81897_338457  0.012884612    0.02576922 6.079152e-09      0.09581288  0.04966505     -0.007588087
# 209239_307893 0.013275289    0.02655058 1.511734e-06      0.71067808 -0.08986896     -0.042189356
# 75590_293847  0.013262380    0.02652476 1.848563e-06      0.27752832  0.03725864     -0.023184643
# 217666_314196 0.008721058    0.01744212 3.880477e-06      0.37038960 -0.05736405     -0.049809589
# 21952_29388   0.028300329    0.05660066 5.600990e-06      0.57027125  0.07304184     -0.078805772
# 69893_298377  0.014361859    0.02872372 6.392233e-06      0.48948751  0.03666609     -0.036445797
# CI_ub_AD CI_ub_Dissection Number_Of_Comparisons Number_of_ Contrasts
# 81897_338457   0.10017180       0.09342541                    31                   31
# 209239_307893 -0.03783078       0.06188699                    38                   38
# 75590_293847   0.08924622       0.08079051                    36                   36
# 217666_314196 -0.02317813       0.01856225                    38                   38
# 21952_29388    0.18397709       0.14306473                    33                   33
# 69893_298377   0.09296354       0.07614911                    38                   38
# Number_of_Coefficients tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity
# 81897_338457                       2               5.617853e-08                  0.0009355331
# 209239_307893                      2               0.000000e+00                  0.0011414147
# 75590_293847                       2               0.000000e+00                  0.0011347241
# 217666_314196                      2               3.020530e-07                  0.0004716520
# 21952_29388                        2               8.557178e-03                  0.0050918059
# 69893_298377                       2               2.185581e-03                  0.0016118909
# QE_CochransQ_Teststat QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity
# 81897_338457               34.39495        0.225096776                    8.633754e-04
# 209239_307893              37.64029        0.394031253                    0.000000e+00
# 75590_293847               28.37195        0.739611913                    0.000000e+00
# 217666_314196              44.83403        0.148344765                    8.535724e-03
# 21952_29388                59.89731        0.001387881                    4.698871e+01
# 69893_298377               54.09401        0.026859516                    3.303253e+01
# H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 81897_338457                            1.000009               2.7738930
# 209239_307893                           1.000000               0.1376006
# 75590_293847                            1.000000               1.1791539
# 217666_314196                           1.000085               0.8023570
# 21952_29388                             1.886391               0.3222289
# 69893_298377                            1.493262               0.4776526
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained       AD_FDR
# 81897_338457                     0.09581288                     99.985791 7.815966e-05
# 209239_307893                    0.71067808                      0.000000 7.922325e-03
# 75590_293847                     0.27752832                    100.000000 7.922325e-03
# 217666_314196                    0.37038960                     99.078290 1.247282e-02
# 21952_29388                      0.57027125                      1.960808 1.369749e-02
# 69893_298377                     0.48948751                      0.000000 1.369749e-02
# Dissection_FDR     QEp_FDR MouseVsRat_EntrezGene.ID Mouse_Symbol
# 81897_338457       0.9372185 0.264491798             81897_338457         Tlr9
# 209239_307893      0.9972600 0.435340708            209239_307893          Gan
# 75590_293847       0.9930209 0.764035865             75590_293847        Dusp9
# 217666_314196      0.9930209 0.182531214            217666_314196       L2hgdh
# 21952_29388        0.9930209 0.002664077              21952_29388        Tnni1
# 69893_298377       0.9930209 0.039643302             69893_298377         Coa7
# Mouse_Genetic.Location Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 81897_338457                Chr9  cM                            Chr9:106099797-106104075(+)
# 209239_307893               Chr8  cM                            Chr8:117884720-117932573(+)
# 75590_293847                ChrX  cM                              ChrX:72683025-72687120(+)
# 217666_314196              Chr12  cM                             Chr12:69737210-69771648(-)
# 21952_29388                 Chr1  cM                            Chr1:135709252-135738727(+)
# 69893_298377                Chr4  cM                            Chr4:108185349-108197915(+)
# Mouse_Name Rat_Symbol Rat_Genetic.Location
# 81897_338457                    toll-like receptor 9       Tlr9                 Chr8
# 209239_307893                giant axonal neuropathy        Gan            Chr19 q12
# 75590_293847          dual specificity phosphatase 9      Dusp9                 ChrX
# 217666_314196     L-2-hydroxyglutarate dehydrogenase     L2hgdh             Chr6 q24
# 21952_29388             troponin I, skeletal, slow 1      Tnni1            Chr13 q13
# 69893_298377  cytochrome c oxidase assembly factor 7       Coa7             Chr5 q34
# Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 81897_338457                                                    NA
# 209239_307893                                                   NA
# 75590_293847                                                    NA
# 217666_314196                                                   NA
# 21952_29388                                                     NA
# 69893_298377                                                    NA
# Rat_Name
# 81897_338457                    toll-like receptor 9
# 209239_307893                              gigaxonin
# 75590_293847          dual specificity phosphatase 9
# 217666_314196     L-2-hydroxyglutarate dehydrogenase
# 21952_29388          troponin I1, slow skeletal type
# 69893_298377  cytochrome c oxidase assembly factor 7


#####################


write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_OrderbyPval.csv")

#Taking a peek

sum(is.na(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl))
#[1] 23
#Very little loss

sum(metaOutputFDR_OrderbyPval$AD_FDR<0.05, na.rm=TRUE)
#[1] 29

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05)]
# [1] "Stx12"    "Psmd1"    "Wdr7"     "Stam"     "L2hgdh"   "Pih1d1"   "Fktn"     "Hk1"     
# [9] "Rpp25l"   "Slc15a2"  "Tlr9"     "Safb2"    "Bcan"     "Dusp9"    "Gan"      "Pnpo"    
# [17] "Stxbp5"   "Coa7"     "Celsr3"   "Tnrc18"   "Kcnip1"   "Slc25a45" "Insc"     "Ido1"    
# [25] "Zcchc7"   "Slco4c1"  "Dmrtb1"   "Tnni1"    "Cdt1"  

hist(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05)], breaks=8)
#All but 3 have abs(Log2FC) <0.10

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05 & abs(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl)>0.10)]
#[1] "Dmrtb1" "Tnni1"  "Cdt1"  

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05 & metaOutputFDR_OrderbyPval$QEp_FDR<0.05)]
#[1] "Bcan"   "Pnpo"   "Stxbp5" "Coa7"   "Tnrc18" "Kcnip1" "Zcchc7" "Dmrtb1" "Tnni1"  "Cdt1"
#So those genes are likely to be highly variable across studies

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05 & metaOutputFDR_OrderbyPval$QEp_CochransQ_pval>0.05)]
# [1] "Stx12"    "Psmd1"    "Wdr7"     "Stam"     "L2hgdh"   "Pih1d1"   "Fktn"     "Hk1"     
# [9] "Rpp25l"   "Slc15a2"  "Tlr9"     "Safb2"    "Dusp9"    "Gan"      "Celsr3"   "Slc25a45"
# [17] "Insc"     "Ido1"     "Slco4c1" 
#Likely to be less variable across studies
#That... totally doesn't seem to match their forest plots.

sum(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05, na.rm=TRUE)
#[1] 1

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05)]
#[1] "Sema4d"

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$AD_FDR<0.05 & metaOutputFDR_OrderbyPval$pval_Dissection<0.05)]
#[1] "Stxbp5"

pdf("Histogram_pvalues_forAD.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_AD)
dev.off()

pdf("Histogram_pvalues_forDissection.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_Dissection)
dev.off()


summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00   24.68   49.87   47.29   70.54   99.39      23

sum(metaOutputFDR_OrderbyPval$QEp_FDR<0.05, na.rm=TRUE)
#[1] 8911

sum(metaOutputFDR_OrderbyPval$QMp_FDR<0.05, na.rm=TRUE)
#[1] 0

plot(metaOutputFDR_OrderbyPval$Log2FC_Dissection_HPCvsPFC~metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl)
#No visible relationship
#Interesting


#Taking a peek at the results for the original DEGs:

OriginalDEGs<-read.csv("OriginalHPC_CTX_DEGs.csv", header=TRUE, stringsAsFactors = FALSE)
str(OriginalDEGs)

OriginalDEGs$Mouse_Symbol
# [1] "Glis1"    "C1qtnf7"  "Acan"     "Tlr9"     "Upk1a"    "Aff3"     "Fbxw10"   "Ccdc160" 
# [9] "Zfp324"   "Cdca7l"   "Fam135b"  "Ccl20"    "Zc3h15"   "Cables1"  "Lrrc1"    "Ccrl2"   
# [17] "Hgsnat"   "Slc25a45" "Iqcg"     "Cdkn2d"   "Mertk"    "Nrsn1"    "Unc13c"   "Tnni1"   
# [25] "Sdhaf3"   "Krt82"    "Stxbp5"   "Hexd"     "Lefty1"   "Clcn2"    "Map3k1"   "Rpp25l"  
# [33] "Dusp9"    "Slc35d2"  "Dcaf10"   "Pea15a"   "Gbp2"     "Or13a28"  "Sim2"     "Relch"   
# [41] "Dmrtb1"   "C2cd2"    "Gpr182"   "Xrcc6"    "Zfp773"   "Ankrd7"   "Zfp41"    "Krt27"   
# [49] "Tbx19"    "Kcnk10"   "Hid1"     "Zfp691"   "Naa50"    "Insrr"    "Fen1"     "Plekha2" 
# [57] "Csad"     "Gpr160"   "Atp6v1b2"

metaOutputFDR_OrderbyPval$Mouse_Symbol[metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol]
# [1] "Rpp25l"   "Nrsn1"    "Relch"    "Tlr9"     "Dusp9"    "Sdhaf3"   "Zfp324"   "Clcn2"   
# [9] "Fen1"     "Stxbp5"   "Naa50"    "Xrcc6"    "Hexd"     "Insrr"    "C2cd2"    "Csad"    
# [17] "Map3k1"   "Ccrl2"    "Slc25a45" "Pea15a"   "Gpr160"   "Hid1"     "Cdca7l"   "Dcaf10"  
# [25] "Atp6v1b2" "Slc35d2"  "Ankrd7"   "Iqcg"     "Cables1"  "Gpr182"   "Lrrc1"    "Zfp691"  
# [33] "Mertk"    "Dmrtb1"   "Tnni1"    "Ccdc160"  "Kcnk10"   "Glis1"    "Unc13c"   "C1qtnf7" 
# [41] "Sim2" 

#Most, but not all, are still in the analysis

write.csv(metaOutputFDR_OrderbyPval[metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol,], "metaOutputFDR_OrderbyPval_forOriginalDEGs.csv")

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol & metaOutputFDR_OrderbyPval$AD_FDR<0.05)]
#[1] "Rpp25l"   "Tlr9"     "Dusp9"    "Stxbp5"   "Slc25a45" "Dmrtb1"   "Tnni1"  

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol & metaOutputFDR_OrderbyPval$pval_AD<0.05)]
# [1] "Rpp25l"   "Nrsn1"    "Relch"    "Tlr9"     "Dusp9"    "Zfp324"   "Clcn2"    "Fen1"    
# [9] "Stxbp5"   "Naa50"    "Hexd"     "Insrr"    "C2cd2"    "Csad"     "Map3k1"   "Slc25a45"
# [17] "Hid1"     "Cdca7l"   "Dcaf10"   "Atp6v1b2" "Slc35d2"  "Ankrd7"   "Iqcg"     "Cables1" 
# [25] "Gpr182"   "Lrrc1"    "Mertk"    "Dmrtb1"   "Tnni1"    "Kcnk10"   "Unc13c"   "C1qtnf7" 


metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol & metaOutputFDR_OrderbyPval$pval_AD>0.05)]
#[1] "Sdhaf3"  "Xrcc6"   "Ccrl2"   "Pea15a"  "Gpr160"  "Zfp691"  "Ccdc160" "Glis1" 


metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Mouse_Symbol%in%OriginalDEGs$Mouse_Symbol & metaOutputFDR_OrderbyPval$pval_Dissection<0.05)]
# [1] "Nrsn1"   "Sdhaf3"  "Zfp324"  "Stxbp5"  "Xrcc6"   "C2cd2"   "Map3k1"  "Ccrl2"   "Pea15a" 
# [10] "Gpr160"  "Cdca7l"  "Slc35d2" "Ankrd7"  "Cables1" "Zfp691"  "Mertk"   "Ccdc160" "Kcnk10" 
# [19] "Glis1"   "Unc13c"


library(metafor)

source("Function_MakeForestPlots_Adjusted.R")

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Rpp25l"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="69961", species="Mouse") 
#convincing

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Tlr9"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="81897", species="Mouse") 
#Noisy but real

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Dusp9"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="75590", species="Mouse") 
#Noisy but maybe convincing?

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Stxbp5"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="78808", species="Mouse") 
#Probably real but small effect 

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Slc25a45"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="107375", species="Mouse") 
#Probably real - a bit noisy, but effect visible

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Dmrtb1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="56296", species="Mouse") 
#Highly variable but probably real

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Tnni1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="21952", species="Mouse") 
#Highly variable but definitely present

#Some of the other top genes
# 81897_338457                    toll-like receptor 9       Tlr9                 Chr8
# 209239_307893                giant axonal neuropathy        Gan            Chr19 q12
# 75590_293847          dual specificity phosphatase 9      Dusp9                 ChrX
# 217666_314196     L-2-hydroxyglutarate dehydrogenase     L2hgdh             Chr6 q24
# 21952_29388             troponin I, skeletal, slow 1      Tnni1            Chr13 q13
# 69893_298377  cytochrome c oxidase assembly factor 7       Coa7             Chr5 q34

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Gan"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="209239", species="Mouse") 
#convincing but a bit noisy

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="L2hgdh"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="217666", species="Mouse") 
#There but small effect size

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Coa7"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="69893", species="Mouse") 
#Visible, small effect size

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Cdt1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="67177", species="Mouse") 
#larger effect size but also quite variable

# [1] "Stx12"    "Psmd1"    "Wdr7"     "Stam"     "L2hgdh"   "Pih1d1"   "Fktn"     "Hk1"     
# [9] "Rpp25l"   "Slc15a2"  "Tlr9"     "Safb2"    "Bcan"     "Dusp9"    "Gan"      "Pnpo"    
# [17] "Stxbp5"   "Coa7"     "Celsr3"   "Tnrc18"   "Kcnip1"   "Slc25a45" "Insc"     "Ido1"    
# [25] "Zcchc7"   "Slco4c1"  "Dmrtb1"   "Tnni1"    "Cdt1"  

#"Ido1"
#This gene encodes indoleamine 2,3-dioxygenase (IDO) - a heme enzyme that catalyzes the first and rate-limiting step in tryptophan catabolism to N-formyl-kynurenine

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Ido1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="15930", species="Mouse")
#Super-duper variable

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Stx12"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="100226", species="Mouse")

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Stam"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="20844", species="Mouse")
# These proteins mediate downstream signaling of cytokine receptors
#Convincing, but very small effect size

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Bcan"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="12032", species="Mouse")
#Mostly DG?

#Hk1
#Hexokinases phosphorylate glucose to produce glucose-6-phosphate, the first step in most glucose metabolism pathways. This gene encodes a ubiquitous form of hexokinase which localizes to the outer membrane of mitochondria
metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Hk1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="15275", species="Mouse")
#Convincing. Effect size is small.

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Safb2"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="224902", species="Mouse")

#The enzyme encoded by this gene catalyzes the terminal, rate-limiting step in the synthesis of pyridoxal 5'-phosphate, also known as vitamin B6
#Pnpo
metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Pnpo"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="103711", species="Mouse")
#Small effect size but actually somewhat convincing

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Kcnip1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="70357", species="Mouse")

metaOutputFDR_annotated$Mouse_EntrezGene.ID[metaOutputFDR_annotated$Mouse_Symbol=="Slco4c1"]
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="227394", species="Mouse")
#Convincing, but a bit noisy

save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_HPCvsCTX/Workspace_MetaRegression_HPC_CTX.RData")

