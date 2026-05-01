#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#In the process, we caught an error in the code (inclusion of GSE84185, which is labeled ACG in the Gemma metadata but is actually the blood data from the same paper as GSE84183) and fixed it
#Updated version: March 10, 2026


### This version is for the PFC


######################

#Grabbing input data and setting the working directory:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/R_Code_And_Workspaces/Meta-analysis")

#Eva directed me to this PFC workspace, but I think it is not the final one
# load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/R_Code_And_Workspaces/Meta-analysis/04112025_PFC_Meta_Run_Reanalysis_Included.RData")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_SophiaEspinoza_Antidepressants_FrontalCortex/ROutput_And_Results/New Results_Reanalysis/April.RData")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff")


######################

#This workspace contains an error - GSE84185 is actually the blood data from the same paper as GSE84183
#We eventually decided to cut it out - but I don't see that in the workspace/code
#And the meta-analysis output in our supplemental file also erroneously reports 17 contrasts - so we must have written up the paper using the wrong workspace/code/output
#I don't know where the final workspace/code/output is with GSE84185 removed, so I'm going to fix our input objects here to remove GSE84185
#Since we're fixing the meta-analysis anyway, I'm also going to make the minimum number of contrasts for inclusion in the meta-analysis the same as what we have for the hippocampus (11)

colnames(MetaAnalysis_FoldChanges)

# [1] "Rat_EntrezGene.ID"             "Mouse_EntrezGene.ID"           "MouseVsRat_EntrezGene.ID"     
# [4] "GSE26836_Amitriptyline"        "GSE84183_fluoxetine"           "GSE118670_Fluoxetine"         
# [7] "GSE28644_fluoxetine"           "GSE93041_ketamine"             "GSE81672_ketamine"            
# [10] "GSE81672_imipramine"           "GSE150264_imipramine"          "GSE84185_fluoxetine"          
# [13] "GSE168172_duloxetine"          "GSE168172_sertraline"          "GSE138802_ketamine"           
# [16] "GSE129359_duloxetine"          "GSE45229_quetiapine_low_dose"  "GSE45229_quetiapine_high_dose"
# [19] "GSE230149_TMS"                 "GSE253280_MDMA"  

MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[,-12]

colnames(MetaAnalysis_FoldChanges)
# [1] "Rat_EntrezGene.ID"             "Mouse_EntrezGene.ID"           "MouseVsRat_EntrezGene.ID"     
# [4] "GSE26836_Amitriptyline"        "GSE84183_fluoxetine"           "GSE118670_Fluoxetine"         
# [7] "GSE28644_fluoxetine"           "GSE93041_ketamine"             "GSE81672_ketamine"            
# [10] "GSE81672_imipramine"           "GSE150264_imipramine"          "GSE168172_duloxetine"         
# [13] "GSE168172_sertraline"          "GSE138802_ketamine"            "GSE129359_duloxetine"         
# [16] "GSE45229_quetiapine_low_dose"  "GSE45229_quetiapine_high_dose" "GSE230149_TMS"                
# [19] "GSE253280_MDMA"  

#correct now

colnames(MetaAnalysis_SV)

MetaAnalysis_SV<-MetaAnalysis_SV[,-12]

colnames(MetaAnalysis_SV)
# [1] "Rat_EntrezGene.ID"             "Mouse_EntrezGene.ID"           "MouseVsRat_EntrezGene.ID"     
# [4] "GSE26836_Amitriptyline"        "GSE84183_fluoxetine"           "GSE118670_Fluoxetine"         
# [7] "GSE28644_fluoxetine"           "GSE93041_ketamine"             "GSE81672_ketamine"            
# [10] "GSE81672_imipramine"           "GSE150264_imipramine"          "GSE168172_duloxetine"         
# [13] "GSE168172_sertraline"          "GSE138802_ketamine"            "GSE129359_duloxetine"         
# [16] "GSE45229_quetiapine_low_dose"  "GSE45229_quetiapine_high_dose" "GSE230149_TMS"                
# [19] "GSE253280_MDMA"   

#Correct now

save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")

#For documentation:
write.csv(MetaAnalysis_FoldChanges, "MetaAnalysis_FoldChanges_PFC.csv")
write.csv(MetaAnalysis_FoldChanges_ForMeta, "MetaAnalysis_FoldChanges_ForMeta_PFC.csv")
write.csv(MetaAnalysis_SV, "MetaAnalysis_SV_PFC.csv")
write.csv(MetaAnalysis_SV_ForMeta, "MetaAnalysis_SV_ForMeta_PFC.csv")


######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

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
  #Then it was made larger to incorporate publication bias (3 stats) 
  #...And even larger to incorporate robustness stats (4 stats)
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 21)
  
  influence_dfbs<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[-c(1:3)]))
  influence_cookd<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[-c(1:3)]))
  influence_TF<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), ncol(MetaAnalysis_FoldChanges_ForMeta[-c(1:3)]))
  
  colnames(influence_dfbs)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)]
  colnames(influence_cookd)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)]
  colnames(influence_TF)<-colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)]
  
  row.names(influence_dfbs)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  row.names(influence_cookd)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  row.names(influence_TF)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    print(i)
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      metaOutput[i, 7]<-TempMeta$k #Metafor output: number of studies (contrasts) included in the analysis - sanity check, should be the same as column 6
      metaOutput[i, 8]<-TempMeta$p #Metafor output: number of coefficients in model
      metaOutput[i, 9]<-TempMeta$tau2 #estimated amount of (residual) heterogeneity
      metaOutput[i, 10]<-TempMeta$se.tau2 #SE of the estimated amount of (residual) heterogeneity
      metaOutput[i, 11]<-TempMeta$QE #test statistic of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 12]<-TempMeta$QEp #p-value of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 13]<-TempMeta$I2 #the I 2 statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes can be attributed to heterogeneity among the true effects
      metaOutput[i, 14]<-TempMeta$H2 #the H2 statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability
      
      #Testing for evidence of publication bias (funnel plot asymmetry)
      skip_to_next <- FALSE
      tryCatch(PubBias<-regtest(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      PubBias<-regtest(TempMeta) #The regression test by Egger et al. (1997) 
      
      metaOutput[i, 15]<-PubBias$zval #the value of the Egger test statistic
      metaOutput[i, 16]<-PubBias$pval #the corresponding Egger test p-value 
      metaOutput[i, 17]<-PubBias$dfs #the corresponding Egger test degrees of freedom
      rm(PubBias)
      }
      
      #Testing for robustness (sensitivity of the meta-analysis to the results of individual contrasts)
      skip_to_next <- FALSE
      tryCatch(Robustness<-leave1out(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      Robustness<-leave1out(TempMeta) #In a leave-one-out analysis, the same model is repeatedly fitted, leaving out one study at a time. By doing so, we can assess how much the results are influenced by each individual study
      
      metaOutput[i, 18]<-min(Robustness$estimate, na.rm=TRUE) #minimum estimated effect in leave-one-out analysis
      metaOutput[i, 19]<-max(Robustness$estimate, na.rm=TRUE) #maximum estimated effect in leave-one-out analysis
      metaOutput[i, 20]<-min(Robustness$pval, na.rm=TRUE) #minimum pval in leave-one-out analysis
      metaOutput[i, 21]<-max(Robustness$pval, na.rm=TRUE) #maximum pval in leave-one-out analysis 
      rm(Robustness)
      }
      
      #Testing for influential datapoints:
      skip_to_next <- FALSE
      tryCatch(InfluenceStats<-influence(TempMeta), error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next){}else{
        
      InfluenceStats<-influence(TempMeta) #Pulling out some stats for influential datapoints (contrasts)
      
      influence_dfbs[i,InfluenceStats$ids]<-InfluenceStats$dfbs$intrcpt
  
      influence_cookd[i,InfluenceStats$ids]<-InfluenceStats$inf$cook.d
      
      influence_TF[i,InfluenceStats$ids]<-InfluenceStats$is.infl
      
      rm(InfluenceStats)
      }
      
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons", "Number_of_ Contrasts", "Number_of_Coefficients", "tau2_ResidualHeterogeneity", "SE_tau2_ResidualHeterogeneity", "QE_CochransQ_Teststat", "QEp_CochransQ_pval", "I2_PercentVar_TrueHeterogeneity", "H2_Ratio_EffectHetero_overSamplVar", "PubBias_Egger_Zstat", "PubBias_Egger_pval", "PubBias_Egger_DF", "Leave1Out_Min_Log2FC","Leave1Out_Max_Log2FC", "Leave1Out_Min_Pval","Leave1Out_Max_Pval")
  
  #The row names for our output are the combined mouse-rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]
  influence_dfbs<<-influence_dfbs
  influence_cookd<<-influence_cookd
  influence_TF<<-influence_TF
  
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  return(influence_dfbs)
  return(influence_cookd)
  return(influence_TF)
  
  #... and provide the user with an update about the newly created object:
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

######################


#Example Usage:

NumberOfComparisons=16
CutOffForNAs=6
#I have 16 comparisons
#5 NA is too many was our original code
#MH: Ideally, this probably should match the HPC - so I changed it to 6 NAs max, with 11 comparisons min

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

######################

#Example Output:


str(metaOutput)
# num [1:15583, 1:21] -0.0068 -0.05028 -0.01136 0.0244 -0.00281 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:21] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# 

head(metaOutput)

# Log2FC_estimate         SE      pval       CI_lb      CI_ub Number_Of_Comparisons
# 23825_114087    -0.006797601 0.01234194 0.5817897 -0.03098736 0.01739216                    16
# 18585_191569    -0.050282937 0.05406125 0.3523139 -0.15624104 0.05567517                    15
# 66514_246307    -0.011362509 0.03010752 0.7058781 -0.07037217 0.04764715                    16
# 20480_65041      0.024402883 0.01819127 0.1797708 -0.01125135 0.06005712                    16
# 13726_25437     -0.002809407 0.02198028 0.8982955 -0.04588996 0.04027115                    16
# 16952_25380      0.104529113 0.09556532 0.2740438 -0.08277548 0.29183371                    15
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 23825_114087                   16                      1                0.000000000
# 18585_191569                   15                      1                0.034622852
# 66514_246307                   16                      1                0.009398246
# 20480_65041                    16                      1                0.002190424
# 13726_25437                    16                      1                0.003170320
# 16952_25380                    15                      1                0.078255969
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 23825_114087                  0.0006861533              6.264427       9.749554e-01
# 18585_191569                  0.0164109329             98.513674       9.132359e-15
# 66514_246307                  0.0050766920             47.894933       2.644191e-05
# 20480_65041                   0.0017510669             27.503192       2.489408e-02
# 13726_25437                   0.0026182397             28.621530       1.798692e-02
# 16952_25380                   0.0486414187             44.273614       5.347797e-05
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar PubBias_Egger_Zstat
# 23825_114087                         0.00000                           1.000000        -0.224277068
# 18585_191569                        87.79823                           8.195535         0.939284909
# 66514_246307                        79.17924                           4.802898         0.007966891
# 20480_65041                         48.60820                           1.945836        -0.433221252
# 13726_25437                         49.43049                           1.977476         0.272109623
# 16952_25380                         73.00519                           3.704416        -0.361673909
# PubBias_Egger_pval PubBias_Egger_DF Leave1Out_Min_Log2FC Leave1Out_Max_Log2FC
# 23825_114087          0.8225417               NA         -0.009474074         -0.002006621
# 18585_191569          0.3475845               NA         -0.069023226         -0.012744959
# 66514_246307          0.9936434               NA         -0.022676290          0.010608597
# 20480_65041           0.6648540               NA          0.007195653          0.032225433
# 13726_25437           0.7855377               NA         -0.011418324          0.005550514
# 16952_25380           0.7175957               NA          0.054765112          0.148668538
# Leave1Out_Min_Pval Leave1Out_Max_Pval
# 23825_114087          0.4506641          0.8760513
# 18585_191569          0.1952691          0.7146329
# 66514_246307          0.4522039          0.9469670
# 20480_65041           0.0858462          0.5809216
# 13726_25437           0.5661581          0.9816990
# 16952_25380           0.1201277          0.5449782

tail(metaOutput)

# Log2FC_estimate         SE       pval       CI_lb        CI_ub Number_Of_Comparisons
# 83673_NA    -0.001042998 0.01737201 0.95212455 -0.03509150  0.033005508                    12
# 93673_NA     0.071301867 0.05154496 0.16657451 -0.02972439  0.172328123                    11
# 93739_NA     0.005242300 0.01195665 0.66106573 -0.01819230  0.028676903                    12
# 93887_NA    -0.030616422 0.02023995 0.13036285 -0.07028600  0.009053157                    12
# 97775_NA     0.030957476 0.02119296 0.14408589 -0.01057996  0.072494916                    12
# 99870_NA    -0.078221628 0.02498038 0.00174021 -0.12718228 -0.029260980                    12
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 83673_NA                   12                      1               4.729851e-06
# 93673_NA                   11                      1               5.775706e-03
# 93739_NA                   12                      1               1.220826e-05
# 93887_NA                   12                      1               0.000000e+00
# 97775_NA                   12                      1               0.000000e+00
# 99870_NA                   12                      1               1.092589e-03
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 83673_NA                  0.0021979446             14.341406         0.21467434
# 93673_NA                  0.0099761851             15.998416         0.09967774
# 93739_NA                  0.0006108196              6.391316         0.84602120
# 93887_NA                  0.0021467157              6.523351         0.83626703
# 97775_NA                  0.0024535679              5.736181         0.89037428
# 99870_NA                  0.0025942165             15.504391         0.16055004
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar PubBias_Egger_Zstat
# 83673_NA                      0.05194076                           1.000520         -0.04131045
# 93673_NA                     24.24320735                           1.320014          2.62835275
# 93739_NA                      0.55460484                           1.005577          1.12462809
# 93887_NA                      0.00000000                           1.000000          0.01116656
# 97775_NA                      0.00000000                           1.000000         -0.24323181
# 99870_NA                     15.42112766                           1.182328          1.19660635
# PubBias_Egger_pval PubBias_Egger_DF Leave1Out_Min_Log2FC Leave1Out_Max_Log2FC
# 83673_NA        0.967048400               NA        -0.0056245805          0.003855966
# 93673_NA        0.008579949               NA         0.0142399379          0.118459201
# 93739_NA        0.260746665               NA         0.0002093189          0.021542383
# 93887_NA        0.991090563               NA        -0.0481637033         -0.022645821
# 97775_NA        0.807825822               NA         0.0235604830          0.034605070
# 99870_NA        0.231460026               NA        -0.0957609890         -0.060583847
# Leave1Out_Min_Pval Leave1Out_Max_Pval
# 83673_NA       0.7470933934         0.97238365
# 93673_NA       0.0433114715         0.68403637
# 93739_NA       0.2046787053         0.98652188
# 93887_NA       0.0957746159         0.27411386
# 97775_NA       0.1093217519         0.50643807
# 99870_NA       0.0008481506         0.00964769


write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")

colnames(metaOutput)

write.csv(influence_dfbs, "influence_dfbs.csv")
write.csv(influence_cookd, "influence_cookd.csv" )
write.csv(influence_TF, "influence_TF.csv")

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
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[22]<-"FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the QEp:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,12], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[23]<-"QEp_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  
  #For the Egger pval
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,16], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[24]<-"Egger_FDR"
  
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
  print(sum(metaOutputFDR_annotated[,24]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]))
  
 
}


FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)

# [1] "metaOutputFDR:"
# num [1:15583, 1:24] -0.0068 -0.05028 -0.01136 0.0244 -0.00281 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:24] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 1
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate          SE         pval
# 11966_117596             117596               11966      0.04871522 0.009256685 1.419497e-07
# 26557_29547               29547               26557     -0.08048271 0.018976741 2.224035e-05
# 16768_297596             297596               16768      0.10684491 0.025237295 2.299678e-05
# 15275_25058               25058               15275      0.04066582 0.010019813 4.938024e-05
# 227394_432363            432363              227394     -0.11618422 0.028947089 5.978050e-05
# 22143_500929             500929               22143      0.04172231 0.010399663 6.023441e-05
# CI_lb       CI_ub Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 11966_117596   0.03057245  0.06685799                    16                   16                      1
# 26557_29547   -0.11767644 -0.04328898                    15                   15                      1
# 16768_297596   0.05738072  0.15630910                    15                   15                      1
# 15275_25058    0.02102735  0.06030429                    16                   16                      1
# 227394_432363 -0.17291947 -0.05944897                    14                   14                      1
# 22143_500929   0.02133935  0.06210527                    14                   14                      1
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 11966_117596                0.000000e+00                  0.0003758395              16.04169
# 26557_29547                 0.000000e+00                  0.0016113693              13.46161
# 16768_297596                1.890744e-06                  0.0028468385              17.85922
# 15275_25058                 0.000000e+00                  0.0004066250              11.46787
# 227394_432363               5.744005e-07                  0.0037060975              13.90855
# 22143_500929                0.000000e+00                  0.0004729750              10.60494
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar
# 11966_117596           0.3792863                     0.000000000                           1.000000
# 26557_29547            0.4905413                     0.000000000                           1.000000
# 16768_297596           0.2132679                     0.015739554                           1.000157
# 15275_25058            0.7187780                     0.000000000                           1.000000
# 227394_432363          0.3803157                     0.003850019                           1.000039
# 22143_500929           0.6438748                     0.000000000                           1.000000
# PubBias_Egger_Zstat PubBias_Egger_pval PubBias_Egger_DF Leave1Out_Min_Log2FC
# 11966_117596           -1.1337063          0.2569178               NA           0.04384071
# 26557_29547            -0.4287161          0.6681299               NA          -0.09112955
# 16768_297596           -0.6453592          0.5186945               NA           0.09577661
# 15275_25058            -1.0429272          0.2969820               NA           0.03655320
# 227394_432363          -1.1726924          0.2409192               NA          -0.15494844
# 22143_500929            0.4669214          0.6405561               NA           0.03308265
# Leave1Out_Max_Log2FC Leave1Out_Min_Pval Leave1Out_Max_Pval         FDR   QEp_FDR
# 11966_117596            0.05219320       4.748323e-08       0.0001265629 0.002212002 0.5064626
# 26557_29547            -0.06942375       4.387497e-06       0.0028695781 0.106984945 0.6102590
# 16768_297596            0.11851389       5.949228e-06       0.0200083873 0.106984945 0.3276177
# 15275_25058             0.04417552       1.887013e-05       0.0012571081 0.106984945 0.8043718
# 227394_432363          -0.10685521       5.471650e-05       0.0002522097 0.106984945 0.5074024
# 22143_500929            0.04470365       3.396796e-05       0.0046424510 0.106984945 0.7447889
# Egger_FDR MouseVsRat_EntrezGene.ID Mouse_Symbol Mouse_Genetic.Location
# 11966_117596  0.9818284             11966_117596     Atp6v1b2               Chr8  cM
# 26557_29547   1.0000000              26557_29547       Homer2               Chr7  cM
# 16768_297596  1.0000000             16768_297596         Lag3               Chr6  cM
# 15275_25058   0.9838184              15275_25058          Hk1              Chr10  cM
# 227394_432363 0.9818284            227394_432363      Slco4c1               Chr1  cM
# 22143_500929  1.0000000             22143_500929       Tuba1b              Chr15  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 11966_117596                               Chr8:69541388-69566370(+)
# 26557_29547                                Chr7:81250229-81356673(-)
# 16768_297596                             Chr6:124881324-124888668(-)
# 15275_25058                               Chr10:62104634-62215687(-)
# 227394_432363                              Chr1:96744918-96800027(-)
# 22143_500929                              Chr15:98829310-98832271(-)
# Mouse_Name Rat_Symbol
# 11966_117596             ATPase, H+ transporting, lysosomal V1 subunit B2   Atp6v1b2
# 26557_29547                                   homer scaffolding protein 2     Homer2
# 16768_297596                                 lymphocyte-activation gene 3       Lag3
# 15275_25058                                                  hexokinase 1        Hk1
# 227394_432363 solute carrier organic anion transporter family, member 4C1    Slco4c1
# 22143_500929                                            tubulin, alpha 1B     Tuba1b
# Rat_Genetic.Location Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 11966_117596             Chr16 p14                                                   NA
# 26557_29547               Chr1 q31                                                   NA
# 16768_297596              Chr4 q42                                                   NA
# 15275_25058              Chr20 q11                                                   NA
# 227394_432363             Chr9 q36                                                   NA
# 22143_500929                  Chr7                                                   NA
# Rat_Name
# 11966_117596                         ATPase H+ transporting V1 subunit B2
# 26557_29547                                      homer scaffold protein 2
# 16768_297596                                      lymphocyte activating 3
# 15275_25058                                                  hexokinase 1
# 227394_432363 solute carrier organic anion transporter family, member 4C1
# 22143_500929                                            tubulin, alpha 1B

#################

#Peeking at the results:

colnames(metaOutputFDR_OrderbyPval)

max(metaOutputFDR_OrderbyPval$Number_Of_Comparisons, na.rm=TRUE)
#[1] 16

min(metaOutputFDR_OrderbyPval$Number_Of_Comparisons, na.rm=TRUE)
#[1] 11

length(metaOutputFDR_OrderbyPval$Log2FC_estimate)
#[1] 15583
sum(is.na(metaOutputFDR_OrderbyPval$Log2FC_estimate)==FALSE)
#[1] 15454

sum(metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 1

sum(metaOutputFDR_OrderbyPval$FDR<0.10, na.rm=TRUE)
#[1] 1

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$FDR<0.05)]
#[1] "Atp6v1b2"


#histogram of I2
hist(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity, breaks=100)
#huge spike at 0

summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.000   0.303  38.836  39.203  68.887  99.590     129 

#For comparison with HPC results:

summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity[c(1:58)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0234  9.9152  9.3123 82.1796 

summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity[-c(1:58)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.4152 39.0024 39.3136 68.9891 99.5897     129 


#histogram of QEp_CochransQ_pval
hist(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval, breaks=100)
#Large spike by 0 - so many genes have significant heterogeneity (not surprising)

sum(metaOutputFDR_OrderbyPval$QEp_FDR<0.05, na.rm=TRUE)
#[1] 6516

sum(metaOutputFDR_OrderbyPval$QEp_FDR[metaOutputFDR_OrderbyPval$FDR<0.05]<0.05, na.rm=TRUE)
#[1] 0

sum(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval[metaOutputFDR_OrderbyPval$FDR<0.05]<0.05, na.rm=TRUE)
#[1] 0
#but not our 1 sig gene lol

#For comparison with HPC:
sum(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval[c(1:58)]<0.05, na.rm=TRUE)
#[1] 7
#Seven of the top 58 genes have nominally significant heterogeneity
sum(metaOutputFDR_OrderbyPval$QEp_FDR[c(1:58)]<0.05, na.rm=TRUE)
#[1] 4 
#four of the top 58 genes have heterogeneity that is sig at FDR<0.05

#histogram of Egger_pval
hist(metaOutputFDR_OrderbyPval$PubBias_Egger_pval, breaks=100)
#Very small little spike at 0

sum(metaOutputFDR_OrderbyPval$Egger_FDR<0.05, na.rm=TRUE)
#[1] 7

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Egger_FDR<0.05)]
#[1] "Pim1"   "Gpr149" "Vgll2"  "Hspa1a" "Cartpt" "Sox1"   "Scg2"

sum(metaOutputFDR_OrderbyPval$PubBias_Egger_pval<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 0

sum(metaOutputFDR_OrderbyPval$Egger_FDR<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 0


library(ggplot2)

TempDF<-data.frame(PubBias=metaOutputFDR_OrderbyPval$Egger_FDR<0.05, Top58=c(rep(TRUE,58), rep(FALSE, (15583-58))), I2=metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)

pdf("ViolinPlot_I2_vs_PubBias.pdf", height=5, width=4)
ggplot(TempDF, aes(x=PubBias, y=I2)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2)
dev.off()

pdf("ViolinPlot_I2_vs_Top58.pdf", height=5, width=4)
ggplot(TempDF, aes(x=Top58, y=I2)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2)
dev.off()


sum(metaOutputFDR_OrderbyPval$Leave1Out_Max_Pval>0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 0

metaOutputFDR_OrderbyPval$Leave1Out_Max_Pval[which(metaOutputFDR_OrderbyPval$FDR<0.05)]
#[1] 0.0001265629

max(metaOutputFDR_OrderbyPval$pval[metaOutputFDR_OrderbyPval$FDR<0.05], na.rm=TRUE)
#[1] 1.419497e-07

sum(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<1.419497e-07, na.rm=TRUE)
#[1] 2

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<1.419497e-07 & metaOutputFDR_OrderbyPval$FDR>0.05)]
#[1] "Gapdhs"

#How many weren't just borderline to begin with?
metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<1.419497e-07 & metaOutputFDR_OrderbyPval$FDR>0.10)]
#[1] "Gapdhs"


#... and how many of those aren't driven by a single study/contrast?
metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<1.419497e-07 & metaOutputFDR_OrderbyPval$Leave1Out_Max_Pval<0.05 & metaOutputFDR_OrderbyPval$FDR>0.10)]
#[1] "Gapdhs"

#Interesting.


#Looking at the influence stats:

str(influence_cookd)
head(influence_cookd)

str(influence_dfbs)
head(influence_dfbs)

str(influence_TF)
head(influence_TF)

boxplot(influence_cookd)
boxplot(influence_dfbs)

Influence_TF_Total<-apply(influence_TF, 2, function(y) sum(y, na.rm=TRUE))
Influence_TF_Total
write.csv(Influence_TF_Total, "Influence_TF_Total.csv")
#The two extreme studies are GSE84183_fluoxetine and GSE28644_fluoxetine

Influence_CooksD_MoreThan1<-apply(influence_cookd, 2, function(y) sum(y>1, na.rm=TRUE))
Influence_CooksD_MoreThan1
write.csv(Influence_CooksD_MoreThan1, "Influence_CooksD_MoreThan1.csv")

Influence_Dfbs_MoreThan1<-apply(influence_dfbs, 2, function(y) sum(abs(y)>1, na.rm=TRUE))
Influence_Dfbs_MoreThan1
write.csv(Influence_Dfbs_MoreThan1, "Influence_Dfbs_MoreThan1.csv")


colnames(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)])

Covariates<-data.frame(row.names=colnames(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)]), ADType=rep(0, 16), Dissection=rep(0, 16), Platform=rep(0, 16), DepressionModel=rep(0, 16), SampleSize=rep(0, 16))

#I think I better fill this in outside of R
write.csv(Covariates, "PFC_Covariates.csv")

#coding:
# ADType: 0.5=NonTrad, -0.5=Trad
# Dissection_ACG: 0.5=ACg, -0.5=PFC, Parietal, Cortex, 0=FC (could have ACg in the mix)
# Dissection_PFC: 0.5=PFC, -0.5=ACg, Parietal, Cortex, 0=FC (could have PFC in the mix)
# Platform: 0.5=microarray, -0.5=RNA-Seq
# DepressionModel: 0.5=depression model included, -0.5=no depression model
# Sample size: treatment + control for that contrast

Covariates<-read.csv("PFC_Covariates.csv", header=TRUE, stringsAsFactors = FALSE)
str(Covariates)

Covariates_AndInfluence<-cbind.data.frame(Covariates, Influence_TF_Total, Influence_CooksD_MoreThan1, Influence_Dfbs_MoreThan1)

cor(as.matrix(Covariates_AndInfluence[,-1]), method="spearman")

write.csv(cor(as.matrix(Covariates_AndInfluence[,-1]), method="spearman"), "CorMatrix_Covariates_spearman.csv")

cor(as.matrix(Covariates_AndInfluence[,-1]), method="pearson")

write.csv(cor(as.matrix(Covariates_AndInfluence[,-1]), method="pearson"), "CorMatrix_Covariates_pearson.csv")


str(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)])

str(as.matrix(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)]))


CorMatrix_Log2FC<-cor(as.matrix(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")

write.csv(CorMatrix_Log2FC, "CorMatrix_Log2FC.csv")

heatmap(CorMatrix_Log2FC)
#col = colorRampPalette(c("blue", "white", "red"))(20)

pdf("Heatmap_CorMatrix_Log2FC_Spearman.pdf", width=8, height=8)
heatmap(CorMatrix_Log2FC, margins = c(20, 20))
dev.off()

pdf("Scatterplot_InfluenceTF_vs_SampleSize.pdf", height=5, width=4)
plot(Influence_TF_Total~SampleSize, data=Covariates_AndInfluence, xlab="Sample Size for Contrast", ylab="# of Genes for which Contrast is Deemed Influential (Outlier)")
dev.off()

pdf("Scatterplot_Influence_CooksD_MoreThan1_vs_SampleSize.pdf", height=5, width=4)
plot(Influence_CooksD_MoreThan1~SampleSize, data=Covariates_AndInfluence, xlab="Sample Size for Contrast", ylab="# of Genes for which Contrast has Cook's D>1 (Outlier)")
dev.off()

pdf("Scatterplot_Influence_Dfbs_MoreThan1_vs_SampleSize.pdf", height=5, width=4)
plot(Influence_Dfbs_MoreThan1~SampleSize, data=Covariates_AndInfluence, xlab="Sample Size for Contrast", ylab="# of Genes for which Contrast has |DFBeta|>1 (Outlier)")
dev.off()


################

#Redoing the heatmap of the top genes


gene_top_50 <- rownames(metaOutputFDR_OrderbyPval)[1:50]
names_genes_top_50<-paste("Mm", metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[1:50], metaOutputFDR_OrderbyPval$Mouse_Symbol[1:50], ";", "Rn", metaOutputFDR_OrderbyPval$Rat_EntrezGene.ID[1:50], metaOutputFDR_OrderbyPval$Rat_Symbol[1:50], sep=" ")

#Rows_Interest<-MetaAnalysis_FoldChanges_ForMeta$MouseVsRat_EntrezGene.ID%in%gene_top_50
#Log2FC_Subsetted <- MetaAnalysis_FoldChanges_ForMeta[Rows_Interest, ]
#Not in the original order:
#cbind(Log2FC_Subsetted$MouseVsRat_EntrezGene.ID, gene_top_50)
#row.names(Log2FC_Subsetted_Matrix)<-Log2FC_Subsetted$MouseVsRat_EntrezGene.ID

Log2FC_Subsetted<-join(data.frame(MouseVsRat_EntrezGene.ID=gene_top_50), MetaAnalysis_FoldChanges_ForMeta, by="MouseVsRat_EntrezGene.ID", type="left")

str(Log2FC_Subsetted)

Log2FC_Subsetted_Matrix <- as.matrix(Log2FC_Subsetted[,-c(1:3)])

row.names(Log2FC_Subsetted_Matrix)<-names_genes_top_50

str(Log2FC_Subsetted_Matrix)

library(pheatmap)
library(dichromat)


pdf("Heatmap_TopMetaGenes_50_color_scale_narrow.pdf",
    height = 11, width = 8.5)

pheatmap(Log2FC_Subsetted_Matrix,
         color         = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 8,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()


##############

#Redoing the forest plot for Atp6v1b2:
#Looks like we don't need to re-do it because it didn't include GSE84185 anyway

##############

#Comparing HPC and PFC results:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions")

list.files()

metaOutputFDR_HPC<-read.csv("metaOutputFDR_orderedByPval_wHeterogeneityPubBiasRobustness.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(metaOutputFDR_HPC)[4:27]<-paste("HPC", colnames(metaOutputFDR_HPC)[4:27], sep="_")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff")

metaOutputFDR_PFC<-metaOutputFDR_OrderbyPval

colnames(metaOutputFDR_PFC)[3:26]<-paste("PFC", colnames(metaOutputFDR_PFC)[3:26], sep="_")

library(plyr)

metaOutputFDR_HPC_PFC<-join(metaOutputFDR_HPC, metaOutputFDR_PFC, by="MouseVsRat_EntrezGene.ID", type="full")

str(metaOutputFDR_HPC_PFC)
#'data.frame':	17488 obs. of  60 variables

pdf("Scatterplot_MetaAnalysisLog2FC_HPC_vs_PFC.pdf", height=5, width=4)
plot(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate~metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate, xlab="HPC: Antidepressant Log2FC", ylab="CTX: Antidepressant Log2FC")
TrendLine<-lm(PFC_Log2FC_estimate~HPC_Log2FC_estimate, data=metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate)==FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate)==FALSE,])
abline(TrendLine, col=2, lwd=3)
dev.off()
#Slight positive correlation

summary.lm(TrendLine)

# Call:
#   lm(formula = PFC_Log2FC_estimate ~ HPC_Log2FC_estimate, data = metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate) == 
#                                                                                          FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate) == 
#                                                                                          FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.36828 -0.01725 -0.00006  0.01740  0.37711 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.0036720  0.0003074   11.95   <2e-16 ***
#   HPC_Log2FC_estimate 0.2856445  0.0064575   44.23   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03691 on 14434 degrees of freedom
# Multiple R-squared:  0.1194,	Adjusted R-squared:  0.1193 
# F-statistic:  1957 on 1 and 14434 DF,  p-value: < 2.2e-16


cor.test(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate, metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate and metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate
# S = 3.5589e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2902152


pdf("Scatterplot_MetaAnalysisLog2FC_HPC_vs_PFC_HPC_FDR10.pdf", height=5, width=4)
plot(PFC_Log2FC_estimate~HPC_Log2FC_estimate,data=metaOutputFDR_HPC_PFC[metaOutputFDR_HPC_PFC$HPC_FDR<0.10,], xlab="HPC: Antidepressant Log2FC", ylab="CTX: Antidepressant Log2FC", pch=16, col="darkgrey")
points(PFC_Log2FC_estimate~HPC_Log2FC_estimate,data=metaOutputFDR_HPC_PFC[metaOutputFDR_HPC_PFC$HPC_FDR<0.10,])
TrendLine<-lm(PFC_Log2FC_estimate~HPC_Log2FC_estimate, data=metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate)==FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate)==FALSE & metaOutputFDR_HPC_PFC$HPC_FDR<0.10,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = PFC_Log2FC_estimate ~ HPC_Log2FC_estimate, data = metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate) == 
#                                                                                          FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate) == 
#                                                                                          FALSE & metaOutputFDR_HPC_PFC$HPC_FDR < 0.1, ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.162829 -0.020470 -0.000796  0.022155  0.114044 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.004396   0.003608   1.218 0.225597    
# HPC_Log2FC_estimate 0.138171   0.036366   3.799 0.000233 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03908 on 116 degrees of freedom
# Multiple R-squared:  0.1107,	Adjusted R-squared:  0.103 
# F-statistic: 14.44 on 1 and 116 DF,  p-value: 0.0002325

cor.test(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR<0.10], metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR<0.10], method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR < 0.1] and metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR < 0.1]
# S = 179286, p-value = 0.0001412
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.345239 


pdf("Scatterplot_MetaAnalysisLog2FC_HPC_vs_PFC_HPC_FDR05.pdf", height=5, width=4)
plot(PFC_Log2FC_estimate~HPC_Log2FC_estimate,data=metaOutputFDR_HPC_PFC[metaOutputFDR_HPC_PFC$HPC_FDR<0.05,], xlab="HPC: Antidepressant Log2FC", ylab="CTX: Antidepressant Log2FC", pch=16, col="darkgrey")
points(PFC_Log2FC_estimate~HPC_Log2FC_estimate,data=metaOutputFDR_HPC_PFC[metaOutputFDR_HPC_PFC$HPC_FDR<0.05,])
TrendLine<-lm(PFC_Log2FC_estimate~HPC_Log2FC_estimate, data=metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate)==FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate)==FALSE & metaOutputFDR_HPC_PFC$HPC_FDR<0.05,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = PFC_Log2FC_estimate ~ HPC_Log2FC_estimate, data = metaOutputFDR_HPC_PFC[is.na(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate) == 
#                                                                                          FALSE & is.na(metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate) == 
#                                                                                          FALSE & metaOutputFDR_HPC_PFC$HPC_FDR < 0.05, ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.162364 -0.024904  0.002058  0.028847  0.114328 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         0.002533   0.006423   0.394   0.6949  
# HPC_Log2FC_estimate 0.127907   0.057740   2.215   0.0313 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04541 on 50 degrees of freedom
# Multiple R-squared:  0.08937,	Adjusted R-squared:  0.07116 
# F-statistic: 4.907 on 1 and 50 DF,  p-value: 0.03133

cor.test(metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR<0.05], metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR<0.05], method="spearman", use="pairwise.complete.obs")
# 
# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_HPC_PFC$PFC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR < 0.05] and metaOutputFDR_HPC_PFC$HPC_Log2FC_estimate[metaOutputFDR_HPC_PFC$HPC_FDR < 0.05]
# S = 17282, p-value = 0.06059
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2622727 

metaOutputFDR_HPC_PFC$Mouse_Symbol[which(metaOutputFDR_HPC_PFC$HPC_FDR<0.05 & metaOutputFDR_HPC_PFC$PFC_pval<0.05)]
#[1] "Tlr9"    "Ccdc160" "Tnni1"   "Stxbp5"  "Dusp9"   "Dmrtb1"  "Insrr"  


##############


save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")


##########


#Re-running fGSEA - code copied from Brain.GMT github and tweaked to fit this analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea", force = TRUE)

library(fgsea)

#This analysis assumes a differential expression (DE) output file structure similar to that produced by the Limma or EdgeR pipelines 
#Rows=all genes included in the DE analysis, columns=gene annotation and DE statistical output
#At least one of the annotation columns must be official gene symbol
#At least one of the columns of differential statistics must include DE effect size (e.g., Log2 Fold Change)

#Read in the full DE results for a condition from the working directory 
#Replace "DEResults.csv" in the code with your file name
DEResults<-metaOutputFDR_OrderbyPval

#Remove rows of DE results that are missing gene symbol annotation or effect size information
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_noNA<-DEResults[is.na(DEResults$Mouse_Symbol)==FALSE & is.na(DEResults$Log2FC_estimate)==FALSE,]

#The analysis only works if there is one effect size (e.g., log2 fold change or Log2FC) per gene symbol.
#One way to deal with multiple effect sizes mapping to the same gene (e.g., multiple transcripts or probes) is to average them:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
DEResults_Log2FC_forGSEA<-tapply(X=DEResults_noNA$Log2FC_estimate, INDEX=DEResults_noNA$Mouse_Symbol, FUN=mean)
names(DEResults_Log2FC_forGSEA)<-names(table(DEResults_noNA$Mouse_Symbol))

#The effect sizes should be ordered from smallest to largest:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_Log2FC_forGSEA_Ranked<-DEResults_Log2FC_forGSEA[order(DEResults_Log2FC_forGSEA)]

str(DEResults_Log2FC_forGSEA_Ranked)
# num [1:14791(1d)] -0.405 -0.369 -0.365 -0.362 -0.334 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:14791] "Cav3" "Slc22a7" "Gpr101" "Slc9a4" ...

#Read in Brain.GMT for your species of interest (this example uses rat)
#If you get a warning about an incomplete line in the .gmt file, just ignore it
setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/R_Code_And_Workspaces/Revisions_Code")
BrainGMT<-gmtPathways("BrainGMTv2_wGO_MouseOrthologs.gmt.txt")


setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff")

#Run fast fGSEA on your ranked, averaged effect sizes:
#This code should be compatible with updated fgsea packages - if you have an updated package, this code will run as fgseaSimple()
GSEA_Results<-fgsea(BrainGMT, DEResults_Log2FC_forGSEA_Ranked, nperm=10000, minSize = 10, maxSize = 1000)

#Pull out the names for the genes that are driving the enrichment of differential expression in each gene set:
GSEA_Results$leadingEdge<-vapply(GSEA_Results$leadingEdge, paste, collapse= ",", character(1L))

#Write out the results:
write.csv(GSEA_Results, "GSEA_Results.csv")

#You can easily view these results in Excel
# Sort by p-value
# padj: false discovery rate (FDR) corrected p-value. This value is normally used to set the threshold for significance (FDR<0.05) 
# ES & NES: Enrichment Score and Normalized Enrichment Score for each gene set. 
# Positive ES & NES values mean that the gene set is enriched with upregulation in response to your variable of interest
# Negative ES & NES values mean that the gene set is enriched with downregulation in response to your variable of interest

# Other aspects of the output can be deciphered by referencing the original GSEA publication: Subramanian et al. 2005
# https://www.pnas.org/doi/10.1073/pnas.0506580102


#Non-directional version:

DEResults_Log2FC_forGSEA_Ranked_NonDirectional<-DEResults_Log2FC_forGSEA[order(abs(DEResults_Log2FC_forGSEA))]

str(DEResults_Log2FC_forGSEA_Ranked_NonDirectional)
# num [1:14791(1d)] -1.74e-06 3.56e-06 4.11e-06 6.08e-06 -1.07e-05 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:14791] "Areg" "Uprt" "Sbds" "Ep300" ...

#Run fast fGSEA on your ranked, averaged effect sizes:
#This code should be compatible with updated fgsea packages - if you have an updated package, this code will run as fgseaSimple()
GSEA_Results_NonDirectional<-fgsea(BrainGMT, DEResults_Log2FC_forGSEA_Ranked_NonDirectional, nperm=10000, minSize = 10, maxSize = 1000)

#Pull out the names for the genes that are driving the enrichment of differential expression in each gene set:
GSEA_Results_NonDirectional$leadingEdge<-vapply(GSEA_Results_NonDirectional$leadingEdge, paste, collapse= ",", character(1L))

#Write out the results:
write.csv(GSEA_Results_NonDirectional, "GSEA_Results_NonDirectional.csv")

