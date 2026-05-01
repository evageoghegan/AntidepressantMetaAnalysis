#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
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
# num [1:16494, 1:21] -0.00287 -0.05409 -0.05024 -0.0274 0.00101 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:21] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)

# Log2FC_estimate         SE       pval       CI_lb        CI_ub Number_Of_Comparisons Number_of_ Contrasts
# 23825_114087   -0.0028695122 0.02317274 0.90144881 -0.04828725  0.042548230                    22                   22
# 18585_191569   -0.0540947362 0.04210825 0.19891150 -0.13662538  0.028435913                    21                   21
# 66514_246307   -0.0502424639 0.02467707 0.04175053 -0.09860863 -0.001876302                    22                   22
# 20480_65041    -0.0273985770 0.01540526 0.07531835 -0.05759233  0.002795179                    22                   22
# 13726_25437     0.0010100471 0.01170981 0.93126252 -0.02194076  0.023960850                    22                   22
# 16952_25380    -0.0001666091 0.02205067 0.99397145 -0.04338513  0.043051913                    21                   21
# Number_of_Coefficients tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 23825_114087                      1               4.976457e-03                  0.0033111895              42.27462
# 18585_191569                      1               2.708820e-02                  0.0115736169              89.02742
# 66514_246307                      1               8.385003e-03                  0.0040210568              74.07272
# 20480_65041                       1               7.315689e-04                  0.0014108748              27.38232
# 13726_25437                       1               2.600137e-06                  0.0006604708              17.20267
# 16952_25380                       1               3.906372e-06                  0.0025304135              22.55717
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar PubBias_Egger_Zstat
# 23825_114087       3.885638e-03                     50.42252543                           2.017045         -0.04746528
# 18585_191569       1.096043e-10                     79.96596923                           4.991507          0.95719031
# 66514_246307       7.703919e-08                     69.62020917                           3.291662         -0.40426289
# 20480_65041        1.585540e-01                     14.67314767                           1.171964          0.18893158
# 13726_25437        6.987495e-01                      0.07946082                           1.000795         -0.39233278
# 16952_25380        3.110576e-01                      0.03678858                           1.000368         -0.98250607
# PubBias_Egger_pval PubBias_Egger_DF Leave1Out_Min_Log2FC Leave1Out_Max_Log2FC Leave1Out_Min_Pval
# 23825_114087          0.9621424               NA         -0.014133363          0.006887216         0.50813591
# 18585_191569          0.3384712               NA         -0.070364797         -0.029516515         0.07990341
# 66514_246307          0.6860194               NA         -0.056494447         -0.015955977         0.02517389
# 20480_65041           0.8501464               NA         -0.031397042         -0.021378385         0.03724838
# 13726_25437           0.6948123               NA         -0.003747755          0.004163859         0.72453216
# 16952_25380           0.3258506               NA         -0.011064885          0.009936567         0.62722762
# Leave1Out_Max_Pval
# 23825_114087          0.9922549
# 18585_191569          0.3058690
# 66514_246307          0.2234328
# 20480_65041           0.1627055
# 13726_25437           0.9859647
# 16952_25380           0.9842227

tail(metaOutput)

# Log2FC_estimate         SE       pval       CI_lb        CI_ub Number_Of_Comparisons Number_of_ Contrasts
# 14934_NA      0.013750373 0.02110330 0.51467576 -0.02761134  0.055112088                    11                   11
# 20568_NA      0.012651714 0.06130017 0.83648665 -0.10749442  0.132797848                    11                   11
# 21991_NA     -0.005753524 0.01652062 0.72764247 -0.03813334  0.026626289                    11                   11
# 246086_NA    -0.030134688 0.01496589 0.04405595 -0.05946729 -0.000802083                    11                   11
# 66166_NA      0.023297459 0.05029703 0.64322280 -0.07528291  0.121877829                    11                   11
# 97848_NA      0.013103586 0.05140867 0.79880761 -0.08765556  0.113862730                    11                   11
# Number_of_Coefficients tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 14934_NA                       1               1.020427e-03                   0.001797956              12.97936
# 20568_NA                       1               2.651182e-02                   0.016856423              35.09146
# 21991_NA                       1               4.188263e-04                   0.001227324              10.69349
# 246086_NA                      1               8.188291e-06                   0.001293316              15.70046
# 66166_NA                       1               1.519579e-02                   0.011221675              41.77991
# 97848_NA                       1               1.911723e-02                   0.011894138              56.47422
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar PubBias_Egger_Zstat
# 14934_NA        2.248280e-01                      24.4666214                           1.323918           0.4827307
# 20568_NA        1.204560e-04                      81.7395037                           5.476302           0.2700050
# 21991_NA        3.818914e-01                      14.1363248                           1.164637          -0.8449282
# 246086_NA       1.085342e-01                       0.1624882                           1.001628           1.5340741
# 66166_NA        8.205878e-06                      69.6979309                           3.300105          -1.4893866
# 97848_NA        1.672967e-08                      80.3845900                           5.098033          -1.0941343
# PubBias_Egger_pval PubBias_Egger_DF Leave1Out_Min_Log2FC Leave1Out_Max_Log2FC Leave1Out_Min_Pval
# 14934_NA           0.6292869               NA          -0.00382560          0.028609953         0.26784097
# 20568_NA           0.7871564               NA          -0.01335833          0.074173678         0.01458794
# 21991_NA           0.3981509               NA          -0.01380291          0.003283356         0.39385587
# 246086_NA          0.1250114               NA          -0.03707846         -0.018123797         0.01621009
# 66166_NA           0.1363856               NA          -0.02108524          0.049259782         0.31287644
# 97848_NA           0.2738961               NA          -0.01815392          0.053788843         0.19566108
# Leave1Out_Max_Pval
# 14934_NA           0.7700130
# 20568_NA           0.9846118
# 21991_NA           0.9972629
# 246086_NA          0.5415932
# 66166_NA           0.8991869
# 97848_NA           0.9967496

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
# num [1:16494, 1:24] -0.00287 -0.05409 -0.05024 -0.0274 0.00101 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:24] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 58
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate         SE         pval       CI_lb       CI_ub
# 230587_298732            298732              230587      0.15627425 0.02334827 2.183672e-11  0.11051248  0.20203601
# 109323_305423            305423              109323     -0.13617634 0.02145280 2.185307e-10 -0.17822305 -0.09412963
# 11595_58968               58968               11595      0.16518348 0.02994748 3.472533e-08  0.10648749  0.22387947
# 81897_338457             338457               81897      0.09635246 0.01757145 4.170683e-08  0.06191306  0.13079186
# 109637_NA                  <NA>              109637      0.09881313 0.01896689 1.890757e-07  0.06163870  0.13598756
# 16764_363220             363220               16764      0.09577155 0.01839426 1.923335e-07  0.05971945  0.13182365
# Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 230587_298732                    14                   14                      1               8.640499e-06
# 109323_305423                    19                   19                      1               4.939213e-06
# 11595_58968                      12                   12                      1               0.000000e+00
# 81897_338457                     16                   16                      1               2.458214e-06
# 109637_NA                        12                   12                      1               0.000000e+00
# 16764_363220                     12                   12                      1               4.231121e-07
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity
# 230587_298732                   0.002280418             18.484915        0.139960137                     0.104959161
# 109323_305423                   0.001947759             39.251105        0.002636425                     0.050458221
# 11595_58968                     0.003836776              8.316288        0.684722267                     0.000000000
# 81897_338457                    0.001284036             17.641701        0.281970520                     0.044823666
# 109637_NA                       0.001461506              4.879375        0.936867183                     0.000000000
# 16764_363220                    0.001454140             17.589711        0.091602094                     0.009804945
# H2_Ratio_EffectHetero_overSamplVar PubBias_Egger_Zstat PubBias_Egger_pval PubBias_Egger_DF
# 230587_298732                           1.001051           1.0533468         0.29218206               NA
# 109323_305423                           1.000505          -0.3283392         0.74265517               NA
# 11595_58968                             1.000000           1.6214231         0.10492693               NA
# 81897_338457                            1.000448          -2.4128073         0.01583019               NA
# 109637_NA                               1.000000          -1.8970603         0.05781999               NA
# 16764_363220                            1.000098           0.7626588         0.44566690               NA
# Leave1Out_Min_Log2FC Leave1Out_Max_Log2FC Leave1Out_Min_Pval Leave1Out_Max_Pval          FDR    QEp_FDR
# 230587_298732           0.14597116            0.1688546       3.021061e-12       1.017630e-07 3.601748e-07 0.19882382
# 109323_305423          -0.14432524           -0.1251177       4.140123e-11       1.243111e-04 1.802223e-06 0.00621217
# 11595_58968             0.15878485            0.1862652       1.060947e-08       2.543959e-06 1.719781e-04 0.73980146
# 81897_338457            0.08429519            0.1019384       1.219218e-08       4.260803e-05 1.719781e-04 0.35712369
# 109637_NA               0.08687121            0.1022430       1.428587e-07       5.869608e-05 5.287249e-04 0.95216509
# 16764_363220            0.08753122            0.1042093       7.640708e-08       2.249258e-04 5.287249e-04 0.13897028
# Egger_FDR MouseVsRat_EntrezGene.ID Mouse_Symbol Mouse_Genetic.Location
# 230587_298732 0.8342545            230587_298732        Glis1               Chr4  cM
# 109323_305423 0.9649520            109323_305423      C1qtnf7               Chr5  cM
# 11595_58968   0.6638530              11595_58968         Acan               Chr7  cM
# 81897_338457  0.3784103             81897_338457         Tlr9               Chr9  cM
# 109637_NA     0.5575106                109637_NA        Upk1a               Chr7  cM
# 16764_363220  0.9006011             16764_363220         Aff3               Chr1  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.                                      Mouse_Name
# 230587_298732                            Chr4:107291788-107492258(+)                       GLIS family zinc finger 1
# 109323_305423                              Chr5:43672881-43776145(+) C1q and tumor necrosis factor related protein 7
# 11595_58968                                Chr7:78703231-78764847(+)                                        aggrecan
# 81897_338457                             Chr9:106099797-106104075(+)                            toll-like receptor 9
# 109637_NA                                  Chr7:30302517-30312159(-)                                    uroplakin 1A
# 16764_363220                               Chr1:38216407-38704036(-)                       AF4/FMR2 family, member 3
# Rat_Symbol Rat_Genetic.Location Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 230587_298732      Glis1             Chr5 q34                                                   NA
# 109323_305423    C1qtnf7            Chr14 q21                                                   NA
# 11595_58968         Acan             Chr1 q31                                                   NA
# 81897_338457        Tlr9                 Chr8                                                   NA
# 109637_NA           <NA>                 <NA>                                                   NA
# 16764_363220        Aff3         Chr9 q21-q22                                                   NA
# Rat_Name
# 230587_298732             GLIS family zinc finger 1
# 109323_305423                 C1q and TNF related 7
# 11595_58968                                aggrecan
# 81897_338457                   toll-like receptor 9
# 109637_NA                                      <NA>
#   16764_363220  ALF transcription elongation factor 3

#################

#Peeking at the results:

colnames(metaOutputFDR_OrderbyPval)

#histogram of I2
hist(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity, breaks=100)

#histogram of I2 for sig genes
hist(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity[metaOutputFDR_OrderbyPval$FDR<0.05], breaks=100)
#Some of the sig genes still have a lot of true heterogeneity (>50%)

sum(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity[metaOutputFDR_OrderbyPval$FDR<0.05]>50, na.rm=TRUE)
#[1] 8


#histogram of QEp_CochransQ_pval
hist(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval, breaks=100)
#Large spike by 0 - so most genes have significant heterogeneity (not surprising)

hist(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval[metaOutputFDR_OrderbyPval$FDR<0.05], breaks=100)

sum(metaOutputFDR_OrderbyPval$QEp_CochransQ_pval[metaOutputFDR_OrderbyPval$FDR<0.05]<0.05, na.rm=TRUE)
#[1] 19

sum(metaOutputFDR_OrderbyPval$QEp_FDR[metaOutputFDR_OrderbyPval$FDR<0.05]<0.05, na.rm=TRUE)
#[1] 14
#So some of the top AD genes have pretty significant heterogeneity

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$QEp_FDR<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05)]
# [1] "C1qtnf7" "Cables1" "Mertk"   "Unc13c"  "Tnni1"   "Stxbp5"  "Rpp25l"  "Pea15a"  "C2cd2"   "Kcnk10"  "Hid1"   
# [12] "Zfp691"  "Plekha2" "Csad" 

#histogram of Egger_pval
hist(metaOutputFDR_OrderbyPval$PubBias_Egger_pval, breaks=100)
#There is a smaller spike by 0 - so maybe some genes??? 
hist(metaOutputFDR_OrderbyPval$PubBias_Egger_pval[metaOutputFDR_OrderbyPval$FDR<0.05], breaks=100)

sum(metaOutputFDR_OrderbyPval$Egger_FDR<0.05, na.rm=TRUE)
#[1] 58 

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Egger_FDR<0.05)]
# [1] "Gdf15"         "Car4"          "Rflnb"         "Gpc1"          "Spon1"         "Pnck"          "Chst7"        
# [8] "Fam43b"        "Clybl"         "Jun"           "Rab26"         "Csgalnact1"    "Calb1"         "Golm1"        
# [15] "Chrnb3"        "Adamts6"       "2610524H06Rik" "Gkn2"          "Podxl2"        "Card6"         "Cd302"      
# [22] "Ccr1l1"        "Gpr83"         "Tdo2"          "Raver1"        "Xpnpep2"       "Pdcd6ip"       "Col27a1"    
# [29] "Tubgcp2"       "Sema3e"        "Fbln2"         "Lhfpl6"        "Ets1"          "Slbp"          "Dock3"      
# [36] "Cdhr2"         "Rps15a"        "Npy"           "Cdk18"         "Adcy10"        "Arl4d"         "Cdh1"       
# [43] "Gphb5"         "Cybb"          "Ptgds"         "Rspo2"         "Rerg"          "Pcdhga2"       "Stra6l"     
# [50] "Egr4"          "Slc22a7"       "Ahdc1"         "Zbtb8b"        "Cdkn2b"        "Tonsl"         "Grm7"       
# [57] "Tnnt1"         "Insl5"   

sum(metaOutputFDR_OrderbyPval$PubBias_Egger_pval<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 6
sum(metaOutputFDR_OrderbyPval$Egger_FDR<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 0
#None of the top AD genes have notable publication bias, but 6 have weak evidence of it (i.e., almost 1 in 10 - so twice what would be expected by random chance)

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$PubBias_Egger_pval<0.05 & metaOutputFDR_OrderbyPval$FDR<0.05)]
#[1] "Tlr9"   "Nrsn1"  "Stxbp5" "Pea15a" "Krt27"  "Hid1" 


#I'm guessing that these genes have greater heterogeneity scores? (is that a given?)

library(ggplot2)

TempDF<-data.frame(PubBias=metaOutputFDR_OrderbyPval$Egger_FDR<0.05, SigAD=metaOutputFDR_OrderbyPval$FDR<0.05, I2=metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)

pdf("ViolinPlot_I2_vs_PubBias.pdf", height=5, width=4)
ggplot(TempDF, aes(x=PubBias, y=I2)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2)
dev.off()

pdf("ViolinPlot_I2_vs_SigAD.pdf", height=5, width=4)
ggplot(TempDF, aes(x=SigAD, y=I2)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2)
dev.off()


sum(metaOutputFDR_OrderbyPval$Leave1Out_Max_Pval>0.05 & metaOutputFDR_OrderbyPval$FDR<0.05, na.rm=TRUE)
#[1] 0
#None of the top genes drop below nominal significance after leaving any of the contrasts out

max(metaOutputFDR_OrderbyPval$pval[metaOutputFDR_OrderbyPval$FDR<0.05], na.rm=TRUE)
#[1] 0.0001753813

sum(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<0.0001753814, na.rm=TRUE)
#[1] 210
#If we consider the same pval threshold that was considered FDR in our meta-analysis
#There are 210 genes that meet that threshold after tossing a contrast out

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<0.0001753814 & metaOutputFDR_OrderbyPval$FDR>0.05)]
# [1] "Hps4"          "Sfxn5"         "Coa7"          "Mstn"          "Hhatl"         "Ints14"        "Cd22"         
# [8] "Gabrr2"        "Or2h1"         "Phtf1"         "C87436"        "Fam228a"       "Ebi3"          "Sema4d"       
# [15] "Abcd4"         "Ifi27l2a"      "Phospho2"      "Tut7"          "Scrn2"         "Ltbp4"         "Knstrn"     
# [22] "Bcan"          "Tfe3"          "Gba2"          "Fscn3"         "Ahnak2"        "L2hgdh"        "Msi2"       
# [29] "Pcgf2"         "Kif12"         "Scart2"        "Commd7"        "Mrps36"        "Aldoart2"      "Zfp14"      
# [36] "1700020L24Rik" "Parvb"         "Eef2"          "Gan"           "Snapc1"        "Arhgap44"      "Dnajc27"    
# [43] "Phex"          "Tesc"          "Mrc1"          "Celsr3"        "Zfp444"        "Ncapg2"        "Lyl1"       
# [50] "Tcp11l1"       "Stx16"         "Twist1"        "Rcc2"          "Ogfod2"        "Acaa2"         "Myh11"      
# [57] "Fbxo28"        "Rfc5"          "Il20"          "Nif3l1"        "Strn"          "Ndfip1"        "Mideas"     
# [64] "Abcc3"         "Foxa3"         "Kcnq1"         "Rwdd2b"        "Krt26"         "Arl9"          "Mospd1"     
# [71] "Osmr"          "Gpt"           "Catsperg1"     "Gdf15"         "Rfx1"          "Nol6"          "Vax2"       
# [78] "Mfsd4b5"       "Tnrc18"        "Hoxb3"         "Pidd1"         "BC034090"      "Trpm4"         "Isg20"      
# [85] "Fbxo45"        "Tmco5"         "Acrbp"         "Cfp"           "Ptpn3"         "Kctd14"        "Tmem150b"   
# [92] "Nsg2"          "Ephx2"         "Gdi2"          "A2ml1"         "Ccl17"         "Myo3b"         "Nudt8"      
# [99] "Ago4"          "Bag4"          "Krt4"          "Garre1"        "St3gal1"       "Dpep1"         "Cimap2"     
# [106] "Fhip2a"        "Celf6"         "Sec22a"        "Fam171a1"      "Castor2"       "Lratd1"        "Cldn14"    
# [113] "Med8"          "Asb14"         "Col11a1"       "Gprin1"        "Hcar1"         "Prodh"         "Sec14l4"   
# [120] "Uso1"          "Chga"          "Lancl1"        "Zfp653"        "Rfc2"          "Klrg2"         "Hyls1"     
# [127] "Slc24a4"       "Melk"          "Tstd1"         "Cadm1"         "Gpihbp1"       "Msantd2"       "Mecr"       
# [134] "Axin2"         "Fmo3"          "Mtmr12"        "Rad9a"         "Cpne1"         "Akap8l"        "Krtap19-5" 
# [141] "Bpifa5"        "Islr2"         "Colgalt2"      "Zfyve28"       "Krt73"         "Slc16a9"       "Ptprd"     
# [148] "3110040N11Rik" "Htra4"         "Shc3"          "Popdc3"        "Tmem232" 

#How many weren't just borderline to begin with?
metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<0.0001753814 & metaOutputFDR_OrderbyPval$FDR>0.10)]
# [1] "Il20"          "Nif3l1"        "Strn"          "Ndfip1"        "Mideas"        "Abcc3"         "Foxa3"        
# [8] "Kcnq1"         "Rwdd2b"        "Krt26"         "Arl9"          "Mospd1"        "Osmr"          "Gpt"          
# [15] "Catsperg1"     "Gdf15"         "Rfx1"          "Nol6"          "Vax2"          "Mfsd4b5"       "Tnrc18"       
# [22] "Hoxb3"         "Pidd1"         "BC034090"      "Trpm4"         "Isg20"         "Fbxo45"        "Tmco5"        
# [29] "Acrbp"         "Cfp"           "Ptpn3"         "Kctd14"        "Tmem150b"      "Nsg2"          "Ephx2"        
# [36] "Gdi2"          "A2ml1"         "Ccl17"         "Myo3b"         "Nudt8"         "Ago4"          "Bag4"         
# [43] "Krt4"          "Garre1"        "St3gal1"       "Dpep1"         "Cimap2"        "Fhip2a"        "Celf6"        
# [50] "Sec22a"        "Fam171a1"      "Castor2"       "Lratd1"        "Cldn14"        "Med8"          "Asb14"        
# [57] "Col11a1"       "Gprin1"        "Hcar1"         "Prodh"         "Sec14l4"       "Uso1"          "Chga"         
# [64] "Lancl1"        "Zfp653"        "Rfc2"          "Klrg2"         "Hyls1"         "Slc24a4"       "Melk"         
# [71] "Tstd1"         "Cadm1"         "Gpihbp1"       "Msantd2"       "Mecr"          "Axin2"         "Fmo3"         
# [78] "Mtmr12"        "Rad9a"         "Cpne1"         "Akap8l"        "Krtap19-5"     "Bpifa5"        "Islr2"        
# [85] "Colgalt2"      "Zfyve28"       "Krt73"         "Slc16a9"       "Ptprd"         "3110040N11Rik" "Htra4"        
# [92] "Shc3"          "Popdc3"        "Tmem232" 


#... and how many of those aren't driven by a single study/contrast?
metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$Leave1Out_Min_Pval<0.0001753814 & metaOutputFDR_OrderbyPval$Leave1Out_Max_Pval<0.05 & metaOutputFDR_OrderbyPval$FDR>0.10)]
# [1] "Il20"      "Nif3l1"    "Strn"      "Ndfip1"    "Mideas"    "Abcc3"     "Foxa3"     "Kcnq1"     "Rwdd2b"   
# [10] "Krt26"     "Arl9"      "Mospd1"    "Osmr"      "Gpt"       "Catsperg1" "Gdf15"     "Rfx1"      "Nol6"     
# [19] "Vax2"      "Mfsd4b5"   "Tnrc18"    "Hoxb3"     "Pidd1"     "BC034090"  "Trpm4"     "Isg20"     "Fbxo45"   
# [28] "Tmco5"     "Acrbp"     "Cfp"       "Ptpn3"     "Kctd14"    "Tmem150b"  "Nsg2"      "Ephx2"     "Gdi2"     
# [37] "A2ml1"     "Ccl17"     "Myo3b"     "Nudt8"     "Ago4"      "Bag4"      "Krt4"      "Garre1"    "St3gal1"  
# [46] "Dpep1"     "Cimap2"    "Celf6"     "Sec22a"    "Fam171a1"  "Castor2"   "Lratd1"    "Cldn14"    "Med8"     
# [55] "Asb14"     "Col11a1"   "Gprin1"    "Hcar1"     "Prodh"     "Sec14l4"   "Uso1"      "Chga"      "Lancl1"   
# [64] "Zfp653"    "Rfc2"      "Melk" 

#Interesting.


#Looking at the influence stats:

str(influence_cookd)
head(influence_cookd)

str(influence_dfbs)
head(influence_dfbs)

str(influence_TF)
head(influence_TF)

#Crap - they are entirely NA values except for the iteration that I ran by hand...
#Time for a lunch break.
#Works now - it just needed to be outputted with a <<-

boxplot(influence_cookd)
boxplot(influence_dfbs)

Influence_TF_Total<-apply(influence_TF, 2, function(y) sum(y, na.rm=TRUE))
write.csv(Influence_TF_Total, "Influence_TF_Total.csv")

Influence_CooksD_MoreThan1<-apply(influence_cookd, 2, function(y) sum(y>1, na.rm=TRUE))
Influence_CooksD_MoreThan1
write.csv(Influence_CooksD_MoreThan1, "Influence_CooksD_MoreThan1.csv")

Influence_Dfbs_MoreThan1<-apply(influence_dfbs, 2, function(y) sum(abs(y)>1, na.rm=TRUE))
Influence_Dfbs_MoreThan1
write.csv(Influence_Dfbs_MoreThan1, "Influence_Dfbs_MoreThan1.csv")


colnames(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)])

Covariates<-data.frame(row.names=colnames(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)]), ADType=rep(0, 22), Dissection=rep(0, 22), Platform=rep(0, 22), DepressionModel=rep(0, 22), SampleSize=rep(0, 22))

#I think I better fill this in outside of R
write.csv(Covariates, "HPC_Covariates.csv")

#coding:
# ADType: 0.5=NonTrad, -0.5=Trad, 0=uncategorized (imipramine+yohimbine)
# Dissection: 0.5=DG, -0.5= Ammon's horn
# Platform: 0.5=microarray, -0.5=RNA-Seq
# DepressionModel: 0.5=depression model included, -0.5=no depression model
# Sample size: treatment + control for that contrast

Covariates<-read.csv("HPC_Covariates.csv", header=TRUE, stringsAsFactors = FALSE)
str(Covariates)

# 'data.frame':	22 obs. of  6 variables:
# $ X              : chr  "GSE123027_ECT" "GSE27532_desipramine" "GSE63469_venlafaxine_high" "GSE63469_venlafaxine_low" ...
# $ ADType         : num  0.5 -0.5 -0.5 -0.5 0.5 0.5 -0.5 -0.5 -0.5 -0.5 ...
# $ Dissection     : num  -0.5 0.5 -0.5 -0.5 -0.5 -0.5 -0.5 0.5 -0.5 0.5 ...
# $ Platform       : num  -0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ DepressionModel: num  0.5 0.5 0.5 0.5 -0.5 0.5 0.5 0.5 -0.5 -0.5 ...
# $ SampleSize     : int  24 16 4 4 24 17 18 32 6 16 ...

Covariates_AndInfluence<-cbind.data.frame(Covariates, Influence_TF_Total, Influence_CooksD_MoreThan1, Influence_Dfbs_MoreThan1)

cor(as.matrix(Covariates_AndInfluence[,-1]), method="spearman")
#                                 ADType  Dissection    Platform DepressionModel  SampleSize Influence_TF_Total
# ADType                      1.00000000 -0.29504722  0.04831022     -0.38263988  0.38380770        -0.19664246
# Dissection                 -0.29504722  1.00000000  0.33062326      0.36991017 -0.02926855         0.29872941
# Platform                    0.04831022  0.33062326  1.00000000     -0.27144836 -0.08386614         0.13570411
# DepressionModel            -0.38263988  0.36991017 -0.27144836      1.00000000 -0.22399454         0.33072031
# SampleSize                  0.38380770 -0.02926855 -0.08386614     -0.22399454  1.00000000        -0.07598594
# Influence_TF_Total         -0.19664246  0.29872941  0.13570411      0.33072031 -0.07598594         1.00000000
# Influence_CooksD_MoreThan1 -0.24886078  0.39389357  0.37622109      0.09239845 -0.22708030         0.78405929
# Influence_Dfbs_MoreThan1   -0.41104000  0.43728847  0.28192659      0.30773369 -0.34997440         0.87376451
# Influence_CooksD_MoreThan1 Influence_Dfbs_MoreThan1
# ADType                                    -0.24886078               -0.4110400
# Dissection                                 0.39389357                0.4372885
# Platform                                   0.37622109                0.2819266
# DepressionModel                            0.09239845                0.3077337
# SampleSize                                -0.22708030               -0.3499744
# Influence_TF_Total                         0.78405929                0.8737645
# Influence_CooksD_MoreThan1                 1.00000000                0.8753184
# Influence_Dfbs_MoreThan1                   0.87531838                1.0000000

write.csv(cor(as.matrix(Covariates_AndInfluence[,-1]), method="spearman"), "CorMatrix_Covariates_spearman.csv")

cor(as.matrix(Covariates_AndInfluence[,-1]), method="pearson")

#                           ADType  Dissection    Platform DepressionModel  SampleSize Influence_TF_Total
# ADType                      1.00000000 -0.28746151  0.04435268      -0.3688007  0.40708889         -0.2187228
# Dissection                 -0.28746151  1.00000000  0.33062326       0.3699102 -0.00851692          0.2762017
# Platform                    0.04435268  0.33062326  1.00000000      -0.2714484 -0.06863737          0.1307803
# DepressionModel            -0.36880069  0.36991017 -0.27144836       1.0000000 -0.24161803          0.2061115
# SampleSize                  0.40708889 -0.00851692 -0.06863737      -0.2416180  1.00000000          0.3752807
# Influence_TF_Total         -0.21872281  0.27620171  0.13078030       0.2061115  0.37528070          1.0000000
# Influence_CooksD_MoreThan1 -0.20162574  0.27867809  0.11958257       0.1592397  0.38893928          0.9638746
# Influence_Dfbs_MoreThan1   -0.23035116  0.28419167  0.12247886       0.1884543  0.35494733          0.9702902
# Influence_CooksD_MoreThan1 Influence_Dfbs_MoreThan1
# ADType                                     -0.2016257               -0.2303512
# Dissection                                  0.2786781                0.2841917
# Platform                                    0.1195826                0.1224789
# DepressionModel                             0.1592397                0.1884543
# SampleSize                                  0.3889393                0.3549473
# Influence_TF_Total                          0.9638746                0.9702902
# Influence_CooksD_MoreThan1                  1.0000000                0.9962154
# Influence_Dfbs_MoreThan1                    0.9962154                1.0000000

write.csv(cor(as.matrix(Covariates_AndInfluence[,-1]), method="pearson"), "CorMatrix_Covariates_pearson.csv")


str(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)])
#'data.frame':	16494 obs. of  22 variables:

str(as.matrix(MetaAnalysis_FoldChanges_ForMeta[,-c(1:3)]))
# num [1:16494, 1:22] -0.0349 0.022 0.0247 -0.0056 0.0258 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "1" "2" "3" "4" ...
# ..$ : chr [1:22] "GSE123027_ECT" "GSE27532_desipramine" "GSE63469_venlafaxine_high" "GSE63469_venlafaxine_low" ...

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


