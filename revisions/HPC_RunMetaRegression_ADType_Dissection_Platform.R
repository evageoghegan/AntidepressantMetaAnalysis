#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include co-variates
#Updated version: March 10, 2026

######################

#Grabbing input data and setting the working directory:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/R_Code_And_Workspaces/Meta-analysis/HPC_November_Final_Meta_Analysis.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions")

######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

ADType<-Covariates$ADType
Dissection<-Covariates$Dissection
Platform<-Covariates$Platform
DepressionModel<-Covariates$DepressionModel

#Platform may make this crash - many of these genes aren't found in the old microarray datasets
  
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
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 32)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    print(i)
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    
    #Making the loop skip to the next gene if there isn't data for one of the co-variates:
    if( (min(table(is.na(effect), Dissection)[1,])<1) | (min(table(is.na(effect), Platform)[1,])<1) ){}else{
    
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~ADType+Dissection+Platform, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~ADType+Dissection+Platform, vi=var)
      metaOutput[i, c(1:4)]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, c(5:8)]<-TempMeta$se #gives standard error
      metaOutput[i, c(9:12)]<-TempMeta$pval #gives pval
      metaOutput[i, c(13:16)]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, c(17:20)]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 21]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      metaOutput[i, 22]<-TempMeta$k #Metafor output: number of studies (contrasts) included in the analysis - sanity check, should be the same as column 6
      metaOutput[i, 23]<-TempMeta$p #Metafor output: number of coefficients in model
      metaOutput[i, 24]<-TempMeta$tau2 #estimated amount of (residual) heterogeneity
      metaOutput[i, 25]<-TempMeta$se.tau2 #SE of the estimated amount of (residual) heterogeneity
      metaOutput[i, 26]<-TempMeta$QE #When moderators are included in the model, this is the QE-test for residual heterogeneity
      metaOutput[i, 27]<-TempMeta$QEp #p-value of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 28]<-TempMeta$I2 #the I 2 statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes can be attributed to heterogeneity among the true effects
      metaOutput[i, 29]<-TempMeta$H2 #the H2 statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability
      metaOutput[i, 30]<-TempMeta$QM #test statistic of the omnibus test of moderators
      metaOutput[i, 31]<-TempMeta$QMp #corresponding p-value of the omnibus test of moderators
      metaOutput[i, 32]<-TempMeta$R2 #amount of heterogeneity accounted for by the moderators
      
      rm(TempMeta)
    }}
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_AD_vs_Ctrl", "Log2FC_ADType_NonTradvsTrad", "Log2FC_Dissection_DG_vs_HPC", "Log2FC_Platform_Microarray_vs_RNASeq",         
                          "SE_AD", "SE_ADType", "SE_Dissection", "SE_Platform",
                          "pval_AD", "pval_ADType","pval_Dissection","pval_Platform",
                          "CI_lb_AD", "CI_lb_ADType","CI_lb_Dissection", "CI_lb_Platform",
                          "CI_ub_AD", "CI_ub_ADType","CI_ub_Dissection", "CI_ub_Platform",
                          "Number_Of_Comparisons", "Number_of_ Contrasts", "Number_of_Coefficients", "tau2_ResidualHeterogeneity", "SE_tau2_ResidualHeterogeneity", "QE_CochransQ_Teststat", "QEp_CochransQ_pval", "I2_PercentVar_TrueHeterogeneity", "H2_Ratio_EffectHetero_overSamplVar", "QM_ModeratorOmnibusTest", "QMp_Pval_ModeratorOmnibusTest", "R2_ModeratorVarianceExplained")
  
  #The row names for our output are the combined mouse-rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]

  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
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

NumberOfComparisons=22
CutOffForNAs=12
#I have 22 comparisons
#12 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

######################

#Example Output:

str(metaOutput)
# num [1:16494, 1:32] -0.05 -0.0366 -0.0803 -0.0416 0.0145 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:32] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "Log2FC_Platform_Microarray_vs_RNASeq" ...

head(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC
# 23825_114087       -0.04999070                -0.093772818                 -0.09736657
# 18585_191569       -0.03658751                 0.048880813                  0.01652209
# 66514_246307       -0.08025966                 0.034401514                 -0.11611055
# 20480_65041        -0.04157152                -0.038625494                 -0.07693336
# 13726_25437         0.01449681                -0.004000575                  0.01674410
# 16952_25380        -0.01208552                 0.011961062                 -0.05915550
# Log2FC_Platform_Microarray_vs_RNASeq      SE_AD  SE_ADType SE_Dissection SE_Platform    pval_AD
# 23825_114087                           0.05737662 0.03008421 0.04688914    0.04180903  0.05980737 0.09657448
# 18585_191569                          -0.03001385 0.07533806 0.10087031    0.10384582  0.14792731 0.62721923
# 66514_246307                           0.05900703 0.03567911 0.04897215    0.05093431  0.06828251 0.02448159
# 20480_65041                            0.01725584 0.02526533 0.03949394    0.04113634  0.05061121 0.09988776
# 13726_25437                           -0.01795030 0.02320585 0.02710234    0.03607448  0.03766115 0.53216468
# 16952_25380                            0.03016861 0.07596878 0.04767810    0.04544115  0.15163088 0.87360166
# pval_ADType pval_Dissection pval_Platform    CI_lb_AD CI_lb_ADType CI_lb_Dissection CI_lb_Platform
# 23825_114087  0.04551282      0.01986750     0.3373789 -0.10895466  -0.18567383       -0.1793108    -0.05984367
# 18585_191569  0.62796672      0.87358839     0.8392164 -0.18424739  -0.14882136       -0.1870120    -0.31994604
# 66514_246307  0.48238554      0.02263060     0.3874998 -0.15018942  -0.06158214       -0.2159400    -0.07482422
# 20480_65041   0.32806909      0.06145547     0.7331419 -0.09109065  -0.11603220       -0.1575591    -0.08194030
# 13726_25437   0.88265060      0.64253783     0.6336281 -0.03098582  -0.05712018       -0.0539606    -0.09176480
# 16952_25380   0.80191371      0.19298312     0.8422934 -0.16098159  -0.08148630       -0.1482185    -0.26702246
# CI_ub_AD CI_ub_ADType CI_ub_Dissection CI_ub_Platform Number_Of_Comparisons Number_of_ Contrasts
# 23825_114087  0.008973267 -0.001871802     -0.015422384      0.1745969                    22                   22
# 18585_191569  0.111072365  0.246582987      0.220056159      0.2599183                    21                   21
# 66514_246307 -0.010329891  0.130385172     -0.016281143      0.1928383                    22                   22
# 20480_65041   0.007947613  0.038781211      0.003692389      0.1164520                    22                   22
# 13726_25437   0.059979443  0.049119029      0.087448786      0.0558642                    22                   22
# 16952_25380   0.136810555  0.105408425      0.029907508      0.3273597                    21                   21
# Number_of_Coefficients tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 23825_114087                      4               1.979661e-03                  0.0023763951              27.10525
# 18585_191569                      4               3.282242e-02                  0.0146805048              84.09091
# 66514_246307                      4               5.733755e-03                  0.0034373543              44.27622
# 20480_65041                       4               1.728943e-03                  0.0021627044              23.47800
# 13726_25437                       4               0.000000e+00                  0.0008069321              16.71250
# 16952_25380                       4               5.177327e-08                  0.0028962005              20.60647
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar
# 23825_114087       7.704242e-02                    2.695567e+01                           1.369032
# 18585_191569       7.134593e-11                    8.149239e+01                           5.403183
# 66514_246307       5.267082e-04                    5.940198e+01                           2.463174
# 20480_65041        1.728835e-01                    2.628356e+01                           1.356549
# 13726_25437        5.429482e-01                    0.000000e+00                           1.000000
# 16952_25380        2.443890e-01                    4.656661e-04                           1.000005
# QM_ModeratorOmnibusTest QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 23825_114087               8.3230314                    0.03978671                      60.21946
# 18585_191569               0.2643844                    0.96658145                       0.00000
# 66514_246307               8.0474013                    0.04504215                      31.61892
# 20480_65041                3.6273301                    0.30461943                       0.00000
# 13726_25437                0.4901768                    0.92104417                     100.00000
# 16952_25380                1.9506801                    0.58270607                      98.67465


tail(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC
# 14934_NA                 NA                          NA                          NA
# 20568_NA                 NA                          NA                          NA
# 21991_NA                 NA                          NA                          NA
# 246086_NA                NA                          NA                          NA
# 66166_NA                 NA                          NA                          NA
# 97848_NA                 NA                          NA                          NA
# Log2FC_Platform_Microarray_vs_RNASeq SE_AD SE_ADType SE_Dissection SE_Platform pval_AD pval_ADType
# 14934_NA                                    NA    NA        NA            NA          NA      NA          NA
# 20568_NA                                    NA    NA        NA            NA          NA      NA          NA
# 21991_NA                                    NA    NA        NA            NA          NA      NA          NA
# 246086_NA                                   NA    NA        NA            NA          NA      NA          NA
# 66166_NA                                    NA    NA        NA            NA          NA      NA          NA
# 97848_NA                                    NA    NA        NA            NA          NA      NA          NA
# pval_Dissection pval_Platform CI_lb_AD CI_lb_ADType CI_lb_Dissection CI_lb_Platform CI_ub_AD CI_ub_ADType
# 14934_NA               NA            NA       NA           NA               NA             NA       NA           NA
# 20568_NA               NA            NA       NA           NA               NA             NA       NA           NA
# 21991_NA               NA            NA       NA           NA               NA             NA       NA           NA
# 246086_NA              NA            NA       NA           NA               NA             NA       NA           NA
# 66166_NA               NA            NA       NA           NA               NA             NA       NA           NA
# 97848_NA               NA            NA       NA           NA               NA             NA       NA           NA
# CI_ub_Dissection CI_ub_Platform Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 14934_NA                NA             NA                    NA                   NA                     NA
# 20568_NA                NA             NA                    NA                   NA                     NA
# 21991_NA                NA             NA                    NA                   NA                     NA
# 246086_NA               NA             NA                    NA                   NA                     NA
# 66166_NA                NA             NA                    NA                   NA                     NA
# 97848_NA                NA             NA                    NA                   NA                     NA
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 14934_NA                          NA                            NA                    NA                 NA
# 20568_NA                          NA                            NA                    NA                 NA
# 21991_NA                          NA                            NA                    NA                 NA
# 246086_NA                         NA                            NA                    NA                 NA
# 66166_NA                          NA                            NA                    NA                 NA
# 97848_NA                          NA                            NA                    NA                 NA
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 14934_NA                               NA                                 NA                      NA
# 20568_NA                               NA                                 NA                      NA
# 21991_NA                               NA                                 NA                      NA
# 246086_NA                              NA                                 NA                      NA
# 66166_NA                               NA                                 NA                      NA
# 97848_NA                               NA                                 NA                      NA
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 14934_NA                             NA                            NA
# 20568_NA                             NA                            NA
# 21991_NA                             NA                            NA
# 246086_NA                            NA                            NA
# 66166_NA                             NA                            NA
# 97848_NA                             NA                            NA


write.csv(metaOutput, "metaOutput_wCovariates_ADType_Dissection_Platform.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wCovariates_ADType_Dissection_Platform.csv")

colnames(metaOutput)

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
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,9], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[33]<-"AD_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the ADType:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,10], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[34]<-"ADType_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  
  #For Dissection:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,11], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[35]<-"Dissection_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  
  #For Platform:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,12], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[36]<-"Platform_FDR"
  
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
  print(sum(metaOutputFDR_annotated[,35]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,35]),]))
  
  
}


FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)

#####################

#Taking a peek

sum(is.na(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl))
#[1] 1413
#Not terrible, but definitely some loss


sum(metaOutputFDR_OrderbyPval$ADType_FDR<0.05, na.rm=TRUE)
#[1] 39

sum(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05, na.rm=TRUE)
#[1] 795

sum(metaOutputFDR_OrderbyPval$Platform_FDR<0.05, na.rm=TRUE)
#[1] 0

#hmmm.... if platform genuinely doesn't seem to matter much, I'm inclined to toss it out
#I am worried about overfitting & there are only 3 RNA-Seq datasets

pdf("Histogram_pvalues_forPlatform.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_Platform)
dev.off()
#Very much not sig

#In contrast:
pdf("Histogram_pvalues_forDissection.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_Dissection)
dev.off()

pdf("Histogram_pvalues_forADType.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_ADType)
dev.off()

pdf("Histogram_pvalues_forAD.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_AD)
dev.off()


