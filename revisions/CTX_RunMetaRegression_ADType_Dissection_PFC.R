#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include co-variates
#Updated version: March 10, 2026

######################

#Grabbing input data and setting the working directory:

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionPFC")

######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

ADType<-Covariates$AD.Type
#Dissection<-Covariates$Dissection_Acg
Dissection<-Covariates$Dissection_PFC
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
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 27)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    print(i)
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    
    #Making the loop skip to the next gene if there isn't data for one of the co-variates:
    if( (min(table(is.na(effect), Dissection)[1,])<1)){}else{
    
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~ADType+Dissection, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~ADType+Dissection, vi=var)
      metaOutput[i, c(1:3)]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, c(4:6)]<-TempMeta$se #gives standard error
      metaOutput[i, c(7:9)]<-TempMeta$pval #gives pval
      metaOutput[i, c(10:12)]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, c(13:15)]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 16]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      metaOutput[i, 17]<-TempMeta$k #Metafor output: number of studies (contrasts) included in the analysis - sanity check, should be the same as column 6
      metaOutput[i, 18]<-TempMeta$p #Metafor output: number of coefficients in model
      metaOutput[i, 19]<-TempMeta$tau2 #estimated amount of (residual) heterogeneity
      metaOutput[i, 20]<-TempMeta$se.tau2 #SE of the estimated amount of (residual) heterogeneity
      metaOutput[i, 21]<-TempMeta$QE #When moderators are included in the model, this is the QE-test for residual heterogeneity
      metaOutput[i, 22]<-TempMeta$QEp #p-value of the test for (residual) heterogeneity (Cochran’s Q-test)
      metaOutput[i, 23]<-TempMeta$I2 #the I 2 statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes can be attributed to heterogeneity among the true effects
      metaOutput[i, 24]<-TempMeta$H2 #the H2 statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability
      metaOutput[i, 25]<-TempMeta$QM #test statistic of the omnibus test of moderators
      metaOutput[i, 26]<-TempMeta$QMp #corresponding p-value of the omnibus test of moderators
      metaOutput[i, 27]<-TempMeta$R2 #amount of heterogeneity accounted for by the moderators
      
      rm(TempMeta)
    }}
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_AD_vs_Ctrl", "Log2FC_ADType_NonTradvsTrad", "Log2FC_Dissection_DG_vs_HPC",         
                          "SE_AD", "SE_ADType", "SE_Dissection", 
                          "pval_AD", "pval_ADType","pval_Dissection",
                          "CI_lb_AD", "CI_lb_ADType","CI_lb_Dissection",
                          "CI_ub_AD", "CI_ub_ADType","CI_ub_Dissection", 
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
# num [1:15583, 1:27] -0.00259 -0.0438 -0.01443 0.0322 -0.0034 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:27] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "SE_AD" ...

head(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC      SE_AD
# 23825_114087      -0.002590123                 0.005267422                  0.02376375 0.01519031
# 18585_191569      -0.043800022                 0.244967273                  0.09274081 0.04875947
# 66514_246307      -0.014429921                 0.015593535                  0.08697802 0.03177824
# 20480_65041        0.032201077                 0.019905217                 -0.05155047 0.01987810
# 13726_25437       -0.003398452                -0.000527824                  0.06538537 0.02415532
# 16952_25380        0.107929606                -0.049617107                  0.16152544 0.09910985
# SE_ADType SE_Dissection   pval_AD pval_ADType pval_Dissection    CI_lb_AD CI_lb_ADType
# 23825_114087 0.03039049    0.02649066 0.8646078  0.86239621       0.3696862 -0.03236258  -0.05429684
# 18585_191569 0.09906462    0.10928903 0.3690323  0.01340581       0.3961134 -0.13936682   0.05080419
# 66514_246307 0.06453890    0.06864354 0.6497699  0.80907887       0.2051206 -0.07671413  -0.11090038
# 20480_65041  0.03923726    0.04069622 0.1052480  0.61194158       0.2052576 -0.00675929  -0.05699840
# 13726_25437  0.04972578    0.05339815 0.8881136  0.99153086       0.2207684 -0.05074201  -0.09798857
# 16952_25380  0.20481654    0.23486314 0.2761584  0.80858532       0.4916147 -0.08632213  -0.45105015
# CI_lb_Dissection   CI_ub_AD CI_ub_ADType CI_ub_Dissection Number_Of_Comparisons
# 23825_114087      -0.02815699 0.02718233   0.06483168       0.07568449                    16
# 18585_191569      -0.12146176 0.05176678   0.43913036       0.30694337                    15
# 66514_246307      -0.04756085 0.04785429   0.14208745       0.22151688                    16
# 20480_65041       -0.13131359 0.07116144   0.09680884       0.02821266                    16
# 13726_25437       -0.03927308 0.04394511   0.09693292       0.17004382                    16
# 16952_25380       -0.29879785 0.30218134   0.35181593       0.62184872                    15
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 23825_114087                   16                      3                0.000000000
# 18585_191569                   15                      3                0.026036327
# 66514_246307                   16                      3                0.009793660
# 20480_65041                    16                      3                0.002259904
# 13726_25437                    16                      3                0.004196055
# 16952_25380                    15                      3                0.086849561
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 23825_114087                   0.001012486              5.375684       9.659655e-01
# 18585_191569                   0.014247745             91.472367       2.559678e-14
# 66514_246307                   0.005745280             43.214441       4.138399e-05
# 20480_65041                    0.002002283             24.062257       3.056219e-02
# 13726_25437                    0.003457749             26.003153       1.698428e-02
# 16952_25380                    0.058471566             35.569909       3.797787e-04
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 23825_114087                         0.00000                           1.000000               0.8887433
# 18585_191569                        80.10496                           5.026379               6.2530580
# 66514_246307                        74.81581                           3.970746               1.6055387
# 20480_65041                         46.58595                           1.872167               1.8905304
# 13726_25437                         50.00724                           2.000289               1.5960695
# 16952_25380                         69.09181                           3.235388               0.6606580
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 23825_114087                     0.6412271                       0.00000
# 18585_191569                     0.0438698                      24.80017
# 66514_246307                     0.4480863                       0.00000
# 20480_65041                      0.3885765                       0.00000
# 13726_25437                      0.4502129                       0.00000
# 16952_25380                      0.7186872                       0.00000

tail(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC      SE_AD
# 83673_NA      -0.003627853               -0.0014281641                 -0.08943444 0.04350502
# 93673_NA       0.123735652                0.0591844434                  0.08200949 0.07683544
# 93739_NA       0.007010170                0.0003397391                  0.02567418 0.01808711
# 93887_NA      -0.023461110                0.0542514845                 -0.05380231 0.03235196
# 97775_NA       0.025793911                0.0319653117                 -0.05488462 0.03183576
# 99870_NA      -0.064512359                0.0305095957                 -0.06913927 0.03338759
# SE_ADType SE_Dissection    pval_AD pval_ADType pval_Dissection    CI_lb_AD CI_lb_ADType
# 83673_NA 0.08711488    0.09345295 0.93354200   0.9869200       0.3385674 -0.08889613  -0.17217019
# 93673_NA 0.15177559    0.17947076 0.10731096   0.6965757       0.6477057 -0.02685904  -0.23829025
# 93739_NA 0.03721784    0.02482371 0.69832821   0.9927167       0.3010145 -0.02843992  -0.07260589
# 93887_NA 0.06294099    0.04839647 0.46833946   0.3887194       0.2662675 -0.08686979  -0.06911058
# 97775_NA 0.06053734    0.05761358 0.41781483   0.5974811       0.3407759 -0.03660304  -0.08668570
# 99870_NA 0.06743838    0.04106448 0.05333163   0.6509758       0.0922444 -0.12995082  -0.10166719
# CI_lb_Dissection     CI_ub_AD CI_ub_ADType CI_ub_Dissection Number_Of_Comparisons
# 83673_NA       -0.2725989 0.0816404258   0.16931386       0.09372998                    12
# 93673_NA       -0.2697467 0.2743303445   0.35665913       0.43376572                    11
# 93739_NA       -0.0229794 0.0424602572   0.07328537       0.07432775                    12
# 93887_NA       -0.1486576 0.0399475725   0.17761355       0.04105302                    12
# 97775_NA       -0.1678052 0.0881908577   0.15061633       0.05803592                    12
# 99870_NA       -0.1496242 0.0009261064   0.16268638       0.01134564                    12
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 83673_NA                   12                      3               7.322458e-03
# 93673_NA                   11                      3               2.241384e-02
# 93739_NA                   12                      3               0.000000e+00
# 93887_NA                   12                      3               0.000000e+00
# 97775_NA                   12                      3               0.000000e+00
# 99870_NA                   12                      3               6.377035e-05
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 83673_NA                   0.010370120             13.683394          0.1340416
# 93673_NA                   0.033957749             12.930462          0.1142647
# 93739_NA                   0.001237680              5.244486          0.8124963
# 93887_NA                   0.005335350              4.780988          0.8529669
# 97775_NA                   0.006830066              4.730213          0.8571637
# 99870_NA                   0.003251804             12.526554          0.1852250
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 83673_NA                      33.0454023                           1.493549               0.9167837
# 93673_NA                      33.5757828                           1.505475               0.3344242
# 93739_NA                       0.0000000                           1.000000               1.1468310
# 93887_NA                       0.0000000                           1.000000               1.7423629
# 97775_NA                       0.0000000                           1.000000               1.0059686
# 99870_NA                       0.7443528                           1.007499               2.8822321
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 83673_NA                     0.6322997                       0.00000
# 93673_NA                     0.8460201                       0.00000
# 93739_NA                     0.5635972                     100.00000
# 93887_NA                     0.4184569                       0.00000
# 97775_NA                     0.6047233                       0.00000
# 99870_NA                     0.2366635                      94.16337

write.csv(metaOutput, "metaOutput_wCovariates_ADType_Dissection_PFC.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wCovariates_ADType_Dissection_PFC.csv")

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
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,7], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[28]<-"AD_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  #For the ADType:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,8], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[29]<-"ADType_FDR"
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
  
  #For Dissection:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,9], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[30]<-"Dissection_FDR"
  
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
  print(sum(metaOutputFDR_annotated[,30]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,35]),]))
  
  
}


FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)

# [1] "metaOutputFDR:"
# num [1:15583, 1:30] -0.00259 -0.0438 -0.01443 0.0322 -0.0034 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:30] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "SE_AD" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 0
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad
# 27993_316317             316317               27993       0.038632675                 0.032375544
# 16008_25662               25662               16008      -0.075796778                -0.011870774
# 12227_29619               29619               12227      -0.019914352                -0.095098626
# 171567_171566            171566              171567      -0.003927398                 0.006973963
# 17901_56781               56781               17901       0.090917433                 0.089712071
# 19165_81751               81751               19165      -0.004605625                 0.036078938
# Log2FC_Dissection_DG_vs_HPC      SE_AD  SE_ADType SE_Dissection    pval_AD pval_ADType
# 27993_316317                   0.01748693 0.01804491 0.03594573    0.03726054 0.03228066   0.3677593
# 16008_25662                   -0.01011976 0.03583582 0.07305765    0.06263254 0.03442038   0.8709239
# 12227_29619                    0.17440930 0.08064074 0.16355568    0.17180388 0.80494576   0.5609406
# 171567_171566                  0.05182454 0.01619665 0.03215828    0.03174131 0.80840655   0.8283146
# 17901_56781                    0.17725989 0.07375021 0.13802215    0.15454335 0.21765965   0.5157031
# 19165_81751                    0.04174675 0.02455064 0.04921018    0.05137093 0.85119261   0.4634608
# pval_Dissection     CI_lb_AD CI_lb_ADType CI_lb_Dissection     CI_ub_AD CI_ub_ADType
# 27993_316317        0.6388445  0.003265296  -0.03807679      -0.05554239  0.074000054   0.10282788
# 16008_25662         0.8716417 -0.146033690  -0.15506113      -0.13287727 -0.005559867   0.13131959
# 12227_29619         0.3100271 -0.177967308  -0.41566187      -0.16232012  0.138138603   0.22546461
# 171567_171566       0.1025288 -0.035672240  -0.05605511      -0.01038730  0.027817444   0.07000303
# 17901_56781         0.2513852 -0.053630328  -0.18080637      -0.12563951  0.235465194   0.36023051
# 19165_81751         0.4164169 -0.052723994  -0.06037124      -0.05893842  0.043512744   0.13252912
# CI_ub_Dissection Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 27993_316317        0.09051626                    15                   15                      3
# 16008_25662         0.11263776                    15                   15                      3
# 12227_29619         0.51113873                    15                   15                      3
# 171567_171566       0.11403637                    16                   16                      3
# 17901_56781         0.48015929                    15                   15                      3
# 19165_81751         0.14243192                    16                   16                      3
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 27993_316317                5.670271e-04                   0.001908872              13.69811
# 16008_25662                 7.663127e-07                   0.005772697              18.64616
# 12227_29619                 5.228524e-02                   0.038122463              26.76517
# 171567_171566               2.277465e-06                   0.001549890              13.25272
# 17901_56781                 2.272342e-02                   0.024109163              33.13409
# 19165_81751                 4.107291e-03                   0.003410511              27.34367
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar
# 27993_316317        0.3204000938                    10.927039458                           1.122675
# 16008_25662         0.0974330451                     0.004235814                           1.000042
# 12227_29619         0.0083514459                    62.787971415                           2.687303
# 171567_171566       0.4284859633                     0.047589394                           1.000476
# 17901_56781         0.0009224534                    57.961427753                           2.378768
# 19165_81751         0.0111594534                    53.705956029                           2.160105
# QM_ModeratorOmnibusTest QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 27993_316317               1.01665871                     0.6014996                      0.000000
# 16008_25662                0.04299857                     0.9787302                      0.000000
# 12227_29619                1.61483771                     0.4460078                      0.000000
# 171567_171566              3.02560535                     0.2202917                     99.195217
# 17901_56781                1.58537478                     0.4526268                      6.198833
# 19165_81751                1.07547548                     0.5840681                      0.000000
# AD_FDR ADType_FDR Dissection_FDR MouseVsRat_EntrezGene.ID Mouse_Symbol
# 27993_316317  0.7315390  0.8803406      0.8946004             27993_316317         Imp4
# 16008_25662   0.7399501  0.9927297      0.9767853              16008_25662       Igfbp2
# 12227_29619   0.9920740  0.9277748      0.7370900              12227_29619         Btg2
# 171567_171566 0.9922124  0.9867927      0.5869010            171567_171566         Nme7
# 17901_56781   0.9039482  0.9150155      0.7027568              17901_56781         Myl1
# 19165_81751   0.9955817  0.9051934      0.7987313              19165_81751        Psen2
# Mouse_Genetic.Location Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 27993_316317                Chr1  cM                              Chr1:34478558-34484828(+)
# 16008_25662                 Chr1  cM                              Chr1:72863662-72891633(+)
# 12227_29619                 Chr1  cM                            Chr1:134002908-134006858(-)
# 171567_171566               Chr1  cM                            Chr1:164135091-164264870(+)
# 17901_56781                 Chr1  cM                              Chr1:66963454-66984563(-)
# 19165_81751                 Chr1  cM                            Chr1:180054569-180091003(-)
# Mouse_Name Rat_Symbol Rat_Genetic.Location
# 27993_316317    IMP4, U3 small nucleolar ribonucleoprotein       Imp4             Chr9 q21
# 16008_25662   insulin-like growth factor binding protein 2     Igfbp2             Chr9 q33
# 12227_29619                BTG anti-proliferation factor 2       Btg2            Chr13 q13
# 171567_171566                     NME/NM23 family member 7       Nme7            Chr13 q23
# 17901_56781                    myosin, light polypeptide 1       Myl1             Chr9 q32
# 19165_81751                                   presenilin 2      Psen2            Chr13 q26
# Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 27993_316317                                                    NA
# 16008_25662                                                     NA
# 12227_29619                                                     NA
# 171567_171566                                                   NA
# 17901_56781                                                     NA
# 19165_81751                                                     NA
# Rat_Name
# 27993_316317    IMP U3 small nucleolar ribonucleoprotein 4
# 16008_25662   insulin-like growth factor binding protein 2
# 12227_29619                BTG anti-proliferation factor 2
# 171567_171566                     NME/NM23 family member 7
# 17901_56781                          myosin, light chain 1
# 19165_81751                                   presenilin 2

#####################

#Taking a peek

sum(is.na(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl))
#[1] 177
#Very little loss


sum(metaOutputFDR_OrderbyPval$ADType_FDR<0.05, na.rm=TRUE)
#[1] 7

sum(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05, na.rm=TRUE)
#[1] 24

sum(metaOutputFDR_OrderbyPval$AD_FDR<0.05, na.rm=TRUE)
#[1] 0


pdf("Histogram_pvalues_forDissection_PFC.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_Dissection)
dev.off()

pdf("Histogram_pvalues_forADType.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_ADType)
dev.off()

pdf("Histogram_pvalues_forAD.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_AD)
dev.off()

#Note that dissection is mislabeled in the data.frame
plot.new()
plot(metaOutputFDR_OrderbyPval$Log2FC_Dissection_DG_vs_HPC ~ metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl)
#Slight negative correlation, not particularly distinct


summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.1596 36.0429 36.9474 64.2922 99.4173     177 

#More heterogeneity remains than when using ACg vs. no ACg as the definition of dissection variation
#But less than without co-variates

############

#Stopped here for now to double check if other definitions of dissection are better


save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionPFC/Workspace_16comparisons6NAcutoff_CorrectComparisons_MetaRegressionADType_PFC.RData")




