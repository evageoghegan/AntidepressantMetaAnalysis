#This code document includes the code for a function that is designed to run a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include co-variates
#Updated version: March 10, 2026

######################

#Grabbing input data and setting the working directory:

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/Workspace_16comparisons6NAcutoff_CorrectComparisons.RData")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionACg")

######################

#Installing and loading relevant code packages:
install.packages("metafor")
library(metafor)
library(plyr) 

######################

ADType<-Covariates$AD.Type
Dissection<-Covariates$Dissection_Acg
#Dissection_PFC<-Covariates$Dissection_PFC
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
# num [1:15583, 1:27] -0.0226 -0.0857 -0.0382 0.0636 -0.0179 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:27] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "SE_AD" ...

head(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC      SE_AD
# 23825_114087       -0.02256256                  0.02248590                 -0.06991655 0.01759828
# 18585_191569       -0.08566260                  0.25588180                 -0.20428971 0.04929313
# 66514_246307       -0.03823857                  0.01611655                 -0.12269161 0.03392172
# 20480_65041         0.06356192                  0.02486790                  0.12593160 0.01591903
# 13726_25437        -0.01786647                 -0.00474856                 -0.06333676 0.02805138
# 16952_25380         0.12831065                 -0.11748200                  0.13147210 0.11269051
# SE_ADType SE_Dissection      pval_AD pval_ADType pval_Dissection    CI_lb_AD CI_lb_ADType
# 23825_114087 0.03090869    0.03476395 1.998116e-01   0.4669231    4.430647e-02 -0.05705456  -0.03809402
# 18585_191569 0.08785041    0.11030486 8.224246e-02   0.0035832    6.401916e-02 -0.18227536   0.08369815
# 66514_246307 0.05991037    0.07324725 2.596329e-01   0.7879211    9.392763e-02 -0.10472392  -0.10130561
# 20480_65041  0.02649443    0.03156994 6.529042e-05   0.3479317    6.636016e-05  0.03236120  -0.02706023
# 13726_25437  0.04996891    0.06265976 5.241773e-01   0.9242908    3.121101e-01 -0.07284617  -0.10268583
# 16952_25380  0.21228349    0.26809540 2.548655e-01   0.5799756    6.238559e-01 -0.09255870  -0.53354999
# CI_lb_Dissection   CI_ub_AD CI_ub_ADType CI_ub_Dissection Number_Of_Comparisons
# 23825_114087      -0.13805264 0.01192945   0.08306582     -0.001780465                    16
# 18585_191569      -0.42048326 0.01095015   0.42806545      0.011903837                    15
# 66514_246307      -0.26625359 0.02824678   0.13353871      0.020870363                    16
# 20480_65041        0.06405566 0.09476265   0.07679603      0.187807539                    16
# 13726_25437       -0.18614763 0.03711323   0.09318871      0.059474112                    16
# 16952_25380       -0.39398523 0.34918000   0.29858599      0.656929434                    15
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 23825_114087                   16                      3                0.000000000
# 18585_191569                   15                      3                0.018882607
# 66514_246307                   16                      3                0.008063469
# 20480_65041                    16                      3                0.000000000
# 13726_25437                    16                      3                0.004414887
# 16952_25380                    15                      3                0.092768323
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 23825_114087                  0.0008473872              2.135559       9.996729e-01
# 18585_191569                  0.0111240539             64.341157       3.607952e-09
# 66514_246307                  0.0049024225             36.073071       5.780459e-04
# 20480_65041                   0.0007474274             10.300332       6.692100e-01
# 13726_25437                   0.0034783672             26.487872       1.460983e-02
# 16952_25380                   0.0611584646             41.851821       3.527228e-05
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 23825_114087                         0.00000                           1.000000               4.1288686
# 18585_191569                        76.31324                           4.221768              10.2362712
# 66514_246307                        75.69462                           4.114316               2.8082290
# 20480_65041                          0.00000                           1.000000              17.2028605
# 13726_25437                         54.29128                           2.187766               1.1132441
# 16952_25380                         70.63696                           3.405642               0.4217519
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 23825_114087                  0.1268900501                       0.00000
# 18585_191569                  0.0059871749                      45.46201
# 66514_246307                  0.2455844212                      14.20241
# 20480_65041                   0.0001838427                     100.00000
# 13726_25437                   0.5731418325                       0.00000
# 16952_25380                   0.8098745121                       0.00000

tail(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC      SE_AD
# 83673_NA        0.02866913                 -0.08121589                0.1455724113 0.03636623
# 93673_NA        0.09294352                 -0.01772635                0.2091084706 0.05073683
# 93739_NA        0.02071452                  0.00263800                0.0328195495 0.02653905
# 93887_NA       -0.01380124                  0.04445625               -0.0001994092 0.03778294
# 97775_NA        0.07236238                 -0.08177500                0.1721551257 0.04106955
# 99870_NA       -0.05344912                  0.01360444                0.0585669578 0.04588678
# SE_ADType SE_Dissection    pval_AD pval_ADType pval_Dissection     CI_lb_AD CI_lb_ADType
# 83673_NA 0.07912470    0.10222812 0.43049498   0.3046894       0.1544477 -0.042607370  -0.23629746
# 93673_NA 0.16065490    0.20637511 0.06697006   0.9121413       0.3109433 -0.006498834  -0.33260416
# 93739_NA 0.03798648    0.05398446 0.43507914   0.9446348       0.5432243 -0.031301068  -0.07181413
# 93887_NA 0.08779914    0.10831060 0.71490464   0.6126178       0.9985310 -0.087854441  -0.12762690
# 97775_NA 0.10097943    0.14091757 0.07807825   0.4180445       0.2218315 -0.008132470  -0.27969104
# 99870_NA 0.08065922    0.09835400 0.24409815   0.8660596       0.5515287 -0.143385557  -0.14448472
# CI_lb_Dissection   CI_ub_AD CI_ub_ADType CI_ub_Dissection Number_Of_Comparisons
# 83673_NA      -0.05479101 0.09994564   0.07386568        0.3459358                    12
# 93673_NA      -0.19537930 0.19238587   0.29715146        0.6135962                    11
# 93739_NA      -0.07298804 0.07273011   0.07709013        0.1386271                    12
# 93887_NA      -0.21248429 0.06025196   0.21653940        0.2120855                    12
# 97775_NA      -0.10403825 0.15285723   0.11614104        0.4483485                    12
# 99870_NA      -0.13420334 0.03648731   0.17169360        0.2513373                    12
# Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 83673_NA                   12                      3               1.590919e-07
# 93673_NA                   11                      3               1.022112e-06
# 93739_NA                   12                      3               1.517429e-05
# 93887_NA                   12                      3               0.000000e+00
# 97775_NA                   12                      3               0.000000e+00
# 99870_NA                   12                      3               2.116545e-03
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 83673_NA                  0.0035394330             12.308137         0.19649011
# 93673_NA                  0.0060696345             12.491141         0.13059870
# 93739_NA                  0.0007033514              5.937605         0.74614645
# 93887_NA                  0.0027512871              6.016860         0.73822972
# 97775_NA                  0.0035936965              4.145239         0.90159564
# 99870_NA                  0.0035988154             15.111421         0.08792055
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 83673_NA                     0.001426996                           1.000014               2.0332578
# 93673_NA                     0.004320569                           1.000043               3.5062691
# 93739_NA                     0.679007596                           1.006836               0.4400885
# 93887_NA                     0.000000000                           1.000000               0.5064911
# 97775_NA                     0.000000000                           1.000000               1.5909419
# 99870_NA                    27.146627105                           1.372620               0.5413268
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 83673_NA                     0.3618126                      96.63643
# 93673_NA                     0.1732301                      99.98230
# 93739_NA                     0.8024833                       0.00000
# 93887_NA                     0.7762773                       0.00000
# 97775_NA                     0.4513686                       0.00000
# 99870_NA                     0.7628732                       0.00000

write.csv(metaOutput, "metaOutput_wCovariates_ADType_Dissection_ACg.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wCovariates_ADType_Dissection_ACg.csv")

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
# num [1:15583, 1:30] -0.0226 -0.0857 -0.0382 0.0636 -0.0179 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15583] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:30] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "SE_AD" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 117
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad
# 27993_316317             316317               27993        0.03268312                  0.03389294
# 16008_25662               25662               16008       -0.04167851                 -0.09393096
# 12227_29619               29619               12227       -0.06330807                 -0.09630279
# 171567_171566            171566              171567       -0.01215976                  0.01486612
# 17901_56781               56781               17901        0.02008330                  0.26317304
# 19165_81751               81751               19165       -0.01419989                  0.03399711
# Log2FC_Dissection_DG_vs_HPC      SE_AD  SE_ADType SE_Dissection   pval_AD  pval_ADType
# 27993_316317                  -0.01744596 0.02417588 0.03522184    0.05289230 0.1764107 0.3359137470
# 16008_25662                    0.16358918 0.04293676 0.09353032    0.11706186 0.3317000 0.3152419539
# 12227_29619                   -0.18565545 0.08715016 0.15839828    0.18677340 0.4675781 0.5432015495
# 171567_171566                 -0.02151420 0.02268070 0.03721274    0.05088904 0.5918700 0.6895321178
# 17901_56781                   -0.20786096 0.03528784 0.07706270    0.05655928 0.5692691 0.0006377015
# 19165_81751                   -0.04097132 0.02657758 0.04890287    0.05610639 0.5931473 0.4869320679
# pval_Dissection    CI_lb_AD CI_lb_ADType CI_lb_Dissection   CI_ub_AD CI_ub_ADType
# 27993_316317     0.7415213922 -0.01470074  -0.03514060      -0.12111295 0.08006698   0.10292649
# 16008_25662      0.1622755233 -0.12583301  -0.27724701      -0.06584786 0.04247599   0.08938509
# 12227_29619      0.3202158471 -0.23411924  -0.40675771      -0.55172459 0.10750309   0.21415213
# 171567_171566    0.6724653624 -0.05661311  -0.05806951      -0.12125488 0.03229359   0.08780175
# 17901_56781      0.0002377573 -0.04907958   0.11213291      -0.31871512 0.08924619   0.41421316
# 19165_81751      0.4652414304 -0.06629099  -0.06185075      -0.15093781 0.03789122   0.12984498
# CI_ub_Dissection Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 27993_316317        0.08622104                    15                   15                      3
# 16008_25662         0.39302622                    15                   15                      3
# 12227_29619         0.18041368                    15                   15                      3
# 171567_171566       0.07822648                    16                   16                      3
# 17901_56781        -0.09700680                    15                   15                      3
# 19165_81751         0.06899518                    16                   16                      3
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat
# 27993_316317                3.311251e-04                   0.001474386              13.92536
# 16008_25662                 2.136394e-06                   0.004115154              16.71912
# 12227_29619                 4.750425e-02                   0.035600988              26.18570
# 171567_171566               6.716784e-04                   0.001674862              15.81887
# 17901_56781                 4.281846e-06                   0.004307945              22.27141
# 19165_81751                 4.072828e-03                   0.003405821              23.94908
# QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar
# 27993_316317          0.30549994                      7.60211754                           1.082276
# 16008_25662           0.16047214                      0.01461036                           1.000146
# 12227_29619           0.01010274                     62.34640604                           2.655789
# 171567_171566         0.25904911                     14.01713011                           1.163022
# 17901_56781           0.03458821                      0.01848405                           1.000185
# 19165_81751           0.03160164                     51.98353737                           2.082619
# QM_ModeratorOmnibusTest QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 27993_316317                0.9431308                  0.6240246542                       0.00000
# 16008_25662                 1.9697901                  0.3734784302                       0.00000
# 12227_29619                 1.5922684                  0.4510693466                       0.00000
# 171567_171566               0.2746218                  0.8716991620                       0.00000
# 17901_56781                17.9426184                  0.0001270018                      99.98232
# 19165_81751                 0.9489512                  0.6222112461                       0.00000
# AD_FDR ADType_FDR Dissection_FDR MouseVsRat_EntrezGene.ID Mouse_Symbol
# 27993_316317  0.5974806  0.8013456     0.90444019             27993_316317         Imp4
# 16008_25662   0.7307905  0.7959195     0.45186762              16008_25662       Igfbp2
# 12227_29619   0.8077377  0.8893932     0.63065109              12227_29619         Btg2
# 171567_171566 0.8671456  0.9402429     0.87079503            171567_171566         Nme7
# 17901_56781   0.8533013  0.1534726     0.01230888              17901_56781         Myl1
# 19165_81751   0.8674872  0.8679940     0.74495039              19165_81751        Psen2
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
#[1] 119
#Very little loss


sum(metaOutputFDR_OrderbyPval$ADType_FDR<0.05, na.rm=TRUE)
#[1] 28

sum(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05, na.rm=TRUE)
#[1] 653

sum(metaOutputFDR_OrderbyPval$AD_FDR<0.05, na.rm=TRUE)
#[1] 117


pdf("Histogram_pvalues_forDissection_ACg.pdf", height=5, width=4)
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
#Strong positive correlation

summary(metaOutputFDR_OrderbyPval$I2_PercentVar_TrueHeterogeneity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.000   0.023  30.747  34.403  61.378  99.412     119


############

#Stopped here for now to double check if other definitions of dissection are better

save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionACg/Workspace_16comparisons6NAcutoff_CorrectComparisons_MetaRegressionADType_ACg.RData")

load("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionACg/Workspace_16comparisons6NAcutoff_CorrectComparisons_MetaRegressionADType_ACg.RData")


#############


#We should probably join that with our already ridiculously massive meta-analysis output to run comparisons

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff")

metaOutputFDR_basic<-read.csv("metaOutputFDR_orderedByPval_wHeterogeneityPubBiasRobustness.csv", header=TRUE, stringsAsFactors = FALSE)

str(metaOutputFDR_basic)
str(metaOutputFDR_OrderbyPval)

metaOutputFDR_cov<-metaOutputFDR_OrderbyPval

colnames(metaOutputFDR_cov)<-paste("wCov", colnames(metaOutputFDR_cov), sep="_")

colnames(metaOutputFDR_cov)[33]<-"MouseVsRat_EntrezGene.ID"  

metaOutputFDR_basic_and_cov<-join(metaOutputFDR_basic, metaOutputFDR_cov, by="MouseVsRat_EntrezGene.ID", type="left")

str(metaOutputFDR_basic_and_cov)
#'data.frame':	15583 obs. of  76 variables:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionACg")

write.csv(metaOutputFDR_basic_and_cov, "metaOutputFDR_basic_and_cov.csv")

pdf("Scatterplot_MetaAnalysisLog2FC_WCov_vs_WO.pdf", height=5, width=4)
plot(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl~metaOutputFDR_basic_and_cov$Log2FC_estimate, xlab="Antidepressant Log2FC: no Covariates", ylab="Antidepressant Log2FC: with Covariates")
abline(a=0, b=1, col="grey", lwd=3)
TrendLine<-lm(wCov_Log2FC_AD_vs_Ctrl ~ Log2FC_estimate, data=metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate)==FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl)==FALSE,])
#I have no idea why R forced me to filter our NA values for that lm.fit - not normal
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_AD_vs_Ctrl ~ Log2FC_estimate, data = metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate) == 
#                                                                                               FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl) == 
#                                                                                               FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37595 -0.01476  0.00035  0.01528  0.24866 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.0006777  0.0002350   2.884  0.00393 ** 
#   Log2FC_estimate 1.2701138  0.0058236 218.097  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02905 on 15363 degrees of freedom
# Multiple R-squared:  0.7559,	Adjusted R-squared:  0.7559 
# F-statistic: 4.757e+04 on 1 and 15363 DF,  p-value: < 2.2e-16


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl, metaOutputFDR_basic_and_cov$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl and metaOutputFDR_basic_and_cov$Log2FC_estimate
# S = 1.0474e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8267469 


#AD Type:

pdf("Scatterplot_MetaAnalysisLog2FC_ADType_vs_Original.pdf", height=5, width=4)
plot(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad~metaOutputFDR_basic_and_cov$Log2FC_estimate, xlab="Antidepressant Log2FC: no Covariates", ylab="Antidepressant Type Log2FC")
#abline(a=0, b=1, col="grey", lwd=3)
TrendLine<-lm(wCov_Log2FC_ADType_NonTradvsTrad~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate)==FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl)==FALSE,])
abline(TrendLine, col=2, lwd=3)
dev.off()
#Slight negative correlation

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_ADType_NonTradvsTrad ~ Log2FC_estimate, 
#      data = metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate) == 
#                                           FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl) == 
#                                           FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.95188 -0.04758 -0.00260  0.04224  1.09571 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.0084702  0.0008447   10.03   <2e-16 ***
#   Log2FC_estimate -0.8490972  0.0209383  -40.55   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1045 on 15363 degrees of freedom
# Multiple R-squared:  0.09669,	Adjusted R-squared:  0.09663 
# F-statistic:  1644 on 1 and 15363 DF,  p-value: < 2.2e-16

cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad, metaOutputFDR_basic_and_cov$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad and metaOutputFDR_basic_and_cov$Log2FC_estimate
# S = 7.3183e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2104996 


plot(metaOutputFDR_basic_and_cov$pval)
#monotonic

head(metaOutputFDR_basic_and_cov$pval)
#[1] 1.419497e-07 2.224035e-05 2.299678e-05 4.938024e-05 5.978050e-05 6.023441e-05
#Yep, still in order by pval

pdf("Scatterplot_MetaAnalysisLog2FC_ADType_vs_Original_Top100.pdf", height=6, width=5)
plot(wCov_Log2FC_ADType_NonTradvsTrad~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),], xlab="Antidepressant Log2FC", ylab="Antidepressant Type Log2FC", pch=16, col="grey")
points(wCov_Log2FC_ADType_NonTradvsTrad~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),])
TrendLine<-lm(wCov_Log2FC_ADType_NonTradvsTrad~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),])
abline(TrendLine, lwd=3, col=2)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_ADType_NonTradvsTrad ~ Log2FC_estimate, 
#      data = metaOutputFDR_basic_and_cov[c(1:100), ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.14755 -0.02234  0.00399  0.03313  0.14289 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     -0.008126   0.005502  -1.477   0.1429  
# Log2FC_estimate -0.172775   0.067193  -2.571   0.0116 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05459 on 98 degrees of freedom
# Multiple R-squared:  0.0632,	Adjusted R-squared:  0.05364 
# F-statistic: 6.612 on 1 and 98 DF,  p-value: 0.01163

cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad[c(1:100)], metaOutputFDR_basic_and_cov$Log2FC_estimate[c(1:100)], method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad[c(1:100)] and metaOutputFDR_basic_and_cov$Log2FC_estimate[c(1:100)]
# S = 191332, p-value = 0.1412
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1481068


pdf("Scatterplot_MetaAnalysisLog2FC_ADType_vs_AD.pdf", height=5, width=4)
plot(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad~metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl, xlab="Antidepressant Log2FC: with Covariates", ylab="Antidepressant Type Log2FC")
#abline(a=0, b=1, col="grey", lwd=3)
TrendLine<-lm(wCov_Log2FC_ADType_NonTradvsTrad~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate)==FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl)==FALSE,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_ADType_NonTradvsTrad ~ wCov_Log2FC_AD_vs_Ctrl, 
#      data = metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate) == 
#                                           FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl) == 
#                                           FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.93536 -0.04847 -0.00138  0.04457  0.97349 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.0091045  0.0008222   11.07   <2e-16 ***
#   wCov_Log2FC_AD_vs_Ctrl -0.7108283  0.0139471  -50.97   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1016 on 15363 degrees of freedom
# Multiple R-squared:  0.1446,	Adjusted R-squared:  0.1446 
# F-statistic:  2598 on 1 and 15363 DF,  p-value: < 2.2e-16


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad, metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad and metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl
# S = 7.3892e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# -0.198906  


pdf("Scatterplot_MetaAnalysisLog2FC_ADType_vs_AD_Top100.pdf", height=6, width=5)
plot(wCov_Log2FC_ADType_NonTradvsTrad~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),], xlab="Antidepressant Log2FC", ylab="Antidepressant Type Log2FC", pch=16, col="grey")
points(wCov_Log2FC_ADType_NonTradvsTrad~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),])
TrendLine<-lm(wCov_Log2FC_ADType_NonTradvsTrad~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),])
abline(TrendLine, lwd=3, col=2)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_ADType_NonTradvsTrad ~ wCov_Log2FC_AD_vs_Ctrl, 
#      data = metaOutputFDR_basic_and_cov[c(1:100), ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.159262 -0.023352  0.004163  0.033434  0.144493 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            -0.008188   0.005568  -1.470   0.1446  
# wCov_Log2FC_AD_vs_Ctrl -0.138373   0.064436  -2.147   0.0342 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05512 on 98 degrees of freedom
# Multiple R-squared:  0.04494,	Adjusted R-squared:  0.0352 
# F-statistic: 4.612 on 1 and 98 DF,  p-value: 0.03422


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad[c(1:100)], metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl[c(1:100)], method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_ADType_NonTradvsTrad[c(1:100)] and metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl[c(1:100)]
# S = 180632, p-value = 0.406
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.08390039 




#Dissection vs. Antidepressant effects:

pdf("Scatterplot_MetaAnalysisLog2FC_Dissection_vs_Original.pdf", height=5, width=4)
plot(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC~metaOutputFDR_basic_and_cov$Log2FC_estimate, xlab="Antidepressant Log2FC: without Covariates", ylab="Dissection Log2FC")
#abline(a=0, b=1, col="grey", lwd=3)
TrendLine<-lm(wCov_Log2FC_Dissection_DG_vs_HPC~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate)==FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl)==FALSE,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_Dissection_DG_vs_HPC ~ Log2FC_estimate, 
#      data = metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate) == 
#                                           FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl) == 
#                                           FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.59342 -0.06312  0.00420  0.06935  1.42780 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -0.006869   0.001210  -5.676 1.41e-08 ***
#   Log2FC_estimate  1.783600   0.029996  59.461  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1496 on 15363 degrees of freedom
# Multiple R-squared:  0.1871,	Adjusted R-squared:  0.187 
# F-statistic:  3536 on 1 and 15363 DF,  p-value: < 2.2e-16

cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC, metaOutputFDR_basic_and_cov$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC and metaOutputFDR_basic_and_cov$Log2FC_estimate
# S = 4.0038e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3377366 


pdf("Scatterplot_MetaAnalysisLog2FC_Dissection_vs_Original_Top100.pdf", height=6, width=5)
plot(wCov_Log2FC_Dissection_DG_vs_HPC~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),], xlab="Antidepressant Log2FC", ylab="Dissection Log2FC", pch=16, col="grey")
points(wCov_Log2FC_Dissection_DG_vs_HPC~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),])
TrendLine<-lm(wCov_Log2FC_Dissection_DG_vs_HPC~Log2FC_estimate, data=metaOutputFDR_basic_and_cov[c(1:100),])
abline(TrendLine, lwd=3, col=2)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_Dissection_DG_vs_HPC ~ Log2FC_estimate, 
#      data = metaOutputFDR_basic_and_cov[c(1:100), ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.304430 -0.040216 -0.007547  0.038308  0.182439 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.009906   0.007333   1.351  0.17984    
# Log2FC_estimate 0.353276   0.089563   3.944  0.00015 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07276 on 98 degrees of freedom
# Multiple R-squared:  0.137,	Adjusted R-squared:  0.1282 
# F-statistic: 15.56 on 1 and 98 DF,  p-value: 0.0001503


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC[c(1:100)], metaOutputFDR_basic_and_cov$Log2FC_estimate[c(1:100)], method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC[c(1:100)] and metaOutputFDR_basic_and_cov$Log2FC_estimate[c(1:100)]
# S = 125186, p-value = 0.01274
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2488089


pdf("Scatterplot_MetaAnalysisLog2FC_Dissection_vs_AD.pdf", height=5, width=4)
plot(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC~metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl, xlab="Antidepressant Log2FC: with Covariates", ylab="Dissection Log2FC")
#abline(a=0, b=1, col="grey", lwd=3)
TrendLine<-lm(wCov_Log2FC_Dissection_DG_vs_HPC~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate)==FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl)==FALSE,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_Dissection_DG_vs_HPC ~ wCov_Log2FC_AD_vs_Ctrl, 
#      data = metaOutputFDR_basic_and_cov[is.na(metaOutputFDR_basic_and_cov$Log2FC_estimate) == 
#                                           FALSE & is.na(metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl) == 
#                                           FALSE, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.36185 -0.04723  0.00344  0.05264  1.06984 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.0105804  0.0009235  -11.46   <2e-16 ***
#   wCov_Log2FC_AD_vs_Ctrl  2.0486863  0.0156659  130.77   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1142 on 15363 degrees of freedom
# Multiple R-squared:  0.5268,	Adjusted R-squared:  0.5267 
# F-statistic: 1.71e+04 on 1 and 15363 DF,  p-value: < 2.2e-16


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC, metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl, method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC and metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl
# S = 2.1479e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6515031 


pdf("Scatterplot_MetaAnalysisLog2FC__Dissection_vs_AD_Top100.pdf", height=6, width=5)
plot(wCov_Log2FC_Dissection_DG_vs_HPC~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),], xlab="Antidepressant Log2FC", ylab="Dissection Log2FC", pch=16, col="grey")
points(wCov_Log2FC_Dissection_DG_vs_HPC~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),])
TrendLine<-lm(wCov_Log2FC_Dissection_DG_vs_HPC~wCov_Log2FC_AD_vs_Ctrl, data=metaOutputFDR_basic_and_cov[c(1:100),])
abline(TrendLine, lwd=3, col=2)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = wCov_Log2FC_Dissection_DG_vs_HPC ~ wCov_Log2FC_AD_vs_Ctrl, 
#      data = metaOutputFDR_basic_and_cov[c(1:100), ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.274806 -0.040200 -0.008931  0.041185  0.171991 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.007799   0.006793   1.148    0.254    
# wCov_Log2FC_AD_vs_Ctrl 0.464795   0.078613   5.912 4.92e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06724 on 98 degrees of freedom
# Multiple R-squared:  0.2629,	Adjusted R-squared:  0.2554 
# F-statistic: 34.96 on 1 and 98 DF,  p-value: 4.921e-08


cor.test(metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC[c(1:100)], metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl[c(1:100)], method="spearman", use="pairwise.complete.obs")

# Spearman's rank correlation rho
# 
# data:  metaOutputFDR_basic_and_cov$wCov_Log2FC_Dissection_DG_vs_HPC[c(1:100)] and metaOutputFDR_basic_and_cov$wCov_Log2FC_AD_vs_Ctrl[c(1:100)]
# S = 90824, p-value = 2.582e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4550015 





#########################

#Checking on top genes:

metaOutputFDR_basic_and_cov$Mouse_Symbol[which(metaOutputFDR_basic_and_cov$FDR<0.05)]
#[1] "Atp6v1b2"


metaOutputFDR_basic_and_cov$Mouse_Symbol[which(metaOutputFDR_basic_and_cov$wCov_AD_FDR<0.05)]

# [1] "Zc3hc1"        "Psmd1"         "Dct"           "Kank4"         "Sprn"          "Gapdhs"       
# [7] "Nudt1"         "Pnpo"          "Ctsz"          "Timm29"        "Snap47"        "Gan"          
# [13] "2310009B15Rik" "Kcnip1"        "Fam163b"       "Pgk1"          "Reep5"         "Stmn2"        
# [19] "Sox30"         "Dnmbp"         "Cacna2d4"      "Fam98a"        "Ank"           "Wdr55"        
# [25] "Paox"          "Snx22"         "Cfl1"          "Tmem177"       "Hnrnpu"        "Atp6v1h"      
# [31] "Gde1"          "Ano4"          "Sobp"          "Fgf1"          "Tango2"        "Yju2b"        
# [37] "Crybb1"        "Necab3"        "Dmrtb1"        "Dnajc6"        "Cdh18"         "Rangap1"      
# [43] "Elf1"          "Serpini1"      "Pmp22"         "Pfkm"          "Gdi2"          "Kcna6"        
# [49] "Mdh1"          "Vamp2"         "Cerox1"        "Nkpd1"         "Ly6e"          "Fam234b"      
# [55] "Chp1"          "Ackr2"         "Polr1a"        "Lrch1"         "Tmx2"          "Ppp1r11"      
# [61] "Atp5f1b"       "Rrp1"          "Nrm"           "Slc48a1"       "Mgst3"         "Asns"         
# [67] "Ckmt1"         "Cltb"          "Dhx8"          "Tubb4b"        "Bcat1"         "Dkk3"         
# [73] "Gmnn"          "Palmd"         "Abhd16b"       "Sostdc1"       "Gnb1"          "Cck"          
# [79] "Etl4"          "Cyth2"         "Alas1"         "Tprkb"         "Rap1gds1"      "Me3"          
# [85] "Mafb"          "Eppk1"         "Rrp7a"         "Fbxo31"        "Gigyf1"        "Mapk8"        
# [91] "Tmem125"       "Fam168a"       "Clpb"          "Micu2"         "Cnot1"         "Anapc1"       
# [97] "Nat14"         "Ndufa6"        "Mtmr7"         "Fchsd2"        "Chrna4"        "Tbca"         
# [103] "Adgrl1"        "Atp2a2"        "Clec18a"       "Psmb5"         "Uba52"         NA             
# [109] "Polr2i"        "Cnbp"          "Stxbp3"        "Snrpg"         "Cox5b"         "Rpl11"        
# [115] "Uap1"          "Ap2s1"         "Ift52"  

#Are the top genes still significant after controlling for covariates?

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_AD_FDR<0.05, na.rm=TRUE)
#[1] 0

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_pval_AD<0.05, na.rm=TRUE)
#[1] 1

metaOutputFDR_basic_and_cov$wCov_pval_AD[which(metaOutputFDR_basic_and_cov$FDR<0.05)]
#[1] 0.01276932
#So Atp6v1b2 is still nominally significant after controlling for covariates


#Is there any evidence that the top genes vary by AD type?

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_ADType_FDR<0.05, na.rm=TRUE)
#[1] 0

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_pval_ADType<0.05, na.rm=TRUE)
#[1] 0

#Is there any evidence that the top genes vary by dissection?

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_Dissection_FDR<0.05, na.rm=TRUE)
#[1] 0

sum(metaOutputFDR_basic_and_cov$FDR<0.05 & metaOutputFDR_basic_and_cov$wCov_pval_Dissection<0.05, na.rm=TRUE)
#[1] 0


pdf("Boxplot_I2_Heterogeneity_Original_vs_MetaRegression.pdf", height=5, width=4)
boxplot(data.frame(Original=metaOutputFDR_basic_and_cov$I2_PercentVar_TrueHeterogeneity, MetaRegression=metaOutputFDR_basic_and_cov$wCov_I2_PercentVar_TrueHeterogeneity), ylab="I2")
dev.off()

pdf("Histogram_ChangeInI2_MetaRegressionVsOriginal.pdf", height=4, width=5)
hist((metaOutputFDR_basic_and_cov$wCov_I2_PercentVar_TrueHeterogeneity-metaOutputFDR_basic_and_cov$I2_PercentVar_TrueHeterogeneity), breaks=100, xlab= "Change in I2: Meta-Regression vs. Original")
dev.off()

summary((metaOutputFDR_basic_and_cov$wCov_I2_PercentVar_TrueHeterogeneity-metaOutputFDR_basic_and_cov$I2_PercentVar_TrueHeterogeneity))
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -92.9750  -9.6273  -1.2283  -4.7399   0.3928  96.6323      218 

save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_EvaGeoghegan_Antidepressants_Hippocampus/ROutput_And_Results/Revisions/PFC/April212025Workspace_16Comparisons_6NAcutoff/MetaRegression_ADType_DissectionACg/Workspace_16comparisons6NAcutoff_CorrectComparisons_MetaRegressionADType_ACg_CompareToOriginal.RData")
