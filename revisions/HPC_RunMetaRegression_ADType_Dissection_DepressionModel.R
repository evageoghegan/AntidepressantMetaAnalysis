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
    if( (min(table(is.na(effect), ADType)[1,c(1,3)])<1) | (min(table(is.na(effect), Dissection)[1,])<1) | (min(table(is.na(effect), DepressionModel)[1,])<1) ){}else{
      
    #Making the loop skip to the next gene if the dissection and depression model are perfectly correlated (redundant predictors)  
     if(round(abs(cor(Dissection[is.na(effect)==FALSE], DepressionModel[is.na(effect)==FALSE])), digits=2)==1){}else{
    
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(yi=effect~ADType+Dissection+DepressionModel, vi=var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(yi=effect~ADType+Dissection+DepressionModel, vi=var)
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
    }}}
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_AD_vs_Ctrl", "Log2FC_ADType_NonTradvsTrad", "Log2FC_Dissection_DG_vs_HPC", "Log2FC_DepressionModel_Yes_vs_No",         
                          "SE_AD", "SE_ADType", "SE_Dissection", "SE_DepressionModel",
                          "pval_AD", "pval_ADType","pval_Dissection","pval_DepressionModel",
                          "CI_lb_AD", "CI_lb_ADType","CI_lb_Dissection", "CI_lb_DepressionModel",
                          "CI_ub_AD", "CI_ub_ADType","CI_ub_Dissection", "CI_ub_DepressionModel",
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
# num [1:16494, 1:32] -0.04776 -0.05277 -0.05996 -0.03335 0.00299 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:32] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "Log2FC_DepressionModel_Yes_vs_No" ...

head(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC
# 23825_114087      -0.047758412               -0.0965841466                -0.092708175
# 18585_191569      -0.052773535                0.0550044818                 0.004958766
# 66514_246307      -0.059958105                0.0439202873                -0.105584581
# 20480_65041       -0.033351287               -0.0409176906                -0.069651012
# 13726_25437        0.002988023                0.0240885552                 0.011518038
# 16952_25380        0.007367163               -0.0007764595                -0.049671322
# Log2FC_DepressionModel_Yes_vs_No      SE_AD  SE_ADType SE_Dissection SE_DepressionModel    pval_AD
# 23825_114087                      0.055133316 0.02570806 0.04323190    0.03677691         0.04865592 0.06320864
# 18585_191569                      0.022589698 0.05195798 0.10933109    0.10068485         0.11611099 0.30977398
# 66514_246307                      0.017428165 0.02554345 0.05278415    0.05032267         0.05537883 0.01891010
# 20480_65041                      -0.008487482 0.02030094 0.03964295    0.03944505         0.04217309 0.10041508
# 13726_25437                       0.047050223 0.01666028 0.03164840    0.03591506         0.02846411 0.85766275
# 16952_25380                      -0.025606877 0.02897366 0.06389624    0.05255493         0.08257836 0.79928619
# pval_ADType pval_Dissection pval_DepressionModel    CI_lb_AD CI_lb_ADType CI_lb_Dissection
# 23825_114087  0.02547685      0.01170800           0.25716106 -0.09814529  -0.18131711      -0.16478960
# 18585_191569  0.61489390      0.96071977           0.84574319 -0.15460931  -0.15928052      -0.19237991
# 66514_246307  0.40536752      0.03589179           0.75298331 -0.11002235  -0.05953475      -0.20421521
# 20480_65041   0.30199923      0.07743389           0.84050036 -0.07314040  -0.11861644      -0.14696188
# 13726_25437   0.44657933      0.74843608           0.09833766 -0.02966552  -0.03794118      -0.05887418
# 16952_25380   0.99030444      0.34459174           0.75649110 -0.04942018  -0.12601078      -0.15267710
# CI_lb_DepressionModel     CI_ub_AD CI_ub_ADType CI_ub_Dissection CI_ub_DepressionModel
# 23825_114087          -0.040230527  0.002628468  -0.01185119     -0.020626748            0.15049716
# 18585_191569          -0.204983655  0.049062238   0.26928948      0.202297438            0.25016305
# 66514_246307          -0.091112348 -0.009893861   0.14737532     -0.006953952            0.12596868
# 20480_65041           -0.091145222  0.006437821   0.03678106      0.007659860            0.07417026
# 13726_25437           -0.008738407  0.035641565   0.08611829      0.081910259            0.10283885
# 16952_25380           -0.187457486  0.064154501   0.12445786      0.053334452            0.13624373
# Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients tau2_ResidualHeterogeneity
# 23825_114087                    22                   22                      4               0.0011115949
# 18585_191569                    21                   21                      4               0.0332058368
# 66514_246307                    22                   22                      4               0.0059357639
# 20480_65041                     22                   22                      4               0.0015450077
# 13726_25437                     22                   22                      4               0.0000000000
# 16952_25380                     21                   21                      4               0.0002277132
# SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval I2_PercentVar_TrueHeterogeneity
# 23825_114087                  0.0019164313              26.55972       8.763190e-02                       17.465613
# 18585_191569                  0.0148379017              85.17623       4.553449e-11                       81.637746
# 66514_246307                  0.0034946808              44.95548       4.205939e-04                       60.386614
# 20480_65041                   0.0020193933              23.46001       1.735263e-01                       24.995658
# 13726_25437                   0.0008104161              14.20737       7.154763e-01                        0.000000
# 16952_25380                   0.0031669386              20.54740       2.471919e-01                        1.893663
# H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest QMp_Pval_ModeratorOmnibusTest
# 23825_114087                           1.211616              10.4915466                    0.01481835
# 18585_191569                           5.445955               0.2582901                    0.96767209
# 66514_246307                           2.524399               7.2096160                    0.06550837
# 20480_65041                            1.333256               3.5984964                    0.30821035
# 13726_25437                            1.000000               2.9953026                    0.39234999
# 16952_25380                            1.019302               1.9342834                    0.58615783
# R2_ModeratorVarianceExplained
# 23825_114087                      77.66292
# 18585_191569                       0.00000
# 66514_246307                      29.20976
# 20480_65041                        0.00000
# 13726_25437                      100.00000
# 16952_25380                        0.00000

tail(metaOutput)

# Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad Log2FC_Dissection_DG_vs_HPC Log2FC_DepressionModel_Yes_vs_No
# 14934_NA        0.080128573                 0.137968319                 0.060585265                      -0.05030491
# 20568_NA        0.001338316                 0.003621083                -0.141097578                       0.09098526
# 21991_NA       -0.018880184                -0.083541904                -0.009056605                      -0.05619302
# 246086_NA      -0.044158921                -0.009591963                -0.012002911                       0.05164779
# 66166_NA        0.050505122                 0.137405971                -0.022416097                       0.16976045
# 97848_NA       -0.007473469                -0.008721970                -0.098712594                       0.11791157
# SE_AD  SE_ADType SE_Dissection SE_DepressionModel   pval_AD pval_ADType pval_Dissection
# 14934_NA  0.06375988 0.14635682    0.06239999         0.07371785 0.2088536   0.3458423       0.3315892
# 20568_NA  0.14021816 0.31029725    0.15768933         0.17762519 0.9923847   0.9906891       0.3709036
# 21991_NA  0.03015246 0.05329583    0.03262789         0.04353691 0.5312117   0.1169948       0.7813405
# 246086_NA 0.06189064 0.13759319    0.07280460         0.06496519 0.4755369   0.9444225       0.8690504
# 66166_NA  0.09844024 0.22351646    0.12156665         0.12621956 0.6079138   0.5387221       0.8537048
# 97848_NA  0.10820971 0.23850051    0.11672636         0.13055919 0.9449381   0.9708278       0.3977339
# pval_DepressionModel    CI_lb_AD CI_lb_ADType CI_lb_Dissection CI_lb_DepressionModel   CI_ub_AD CI_ub_ADType
# 14934_NA             0.4949874 -0.04483849   -0.1488858      -0.06171646           -0.19478925 0.20509563   0.42482241
# 20568_NA             0.6084888 -0.27348422   -0.6045503      -0.45016298           -0.25715371 0.27616086   0.61179251
# 21991_NA             0.1968082 -0.07797791   -0.1879998      -0.07300609           -0.14152379 0.04021755   0.02091599
# 246086_NA            0.4266093 -0.16546234   -0.2792697      -0.15469730           -0.07568164 0.07714450   0.26008574
# 66166_NA             0.1786377 -0.14243420   -0.3006782      -0.26068236           -0.07762534 0.24344445   0.57549018
# 97848_NA             0.3664583 -0.21956060   -0.4761744      -0.32749206           -0.13797974 0.20461366   0.45873044
# CI_ub_Dissection CI_ub_DepressionModel Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 14934_NA        0.18288699            0.09417943                    11                   11                      4
# 20568_NA        0.16796782            0.43912423                    11                   11                      4
# 21991_NA        0.05489288            0.02913775                    11                   11                      4
# 246086_NA       0.13069147            0.17897723                    11                   11                      4
# 66166_NA        0.21585017            0.41714625                    11                   11                      4
# 97848_NA        0.13006687            0.37380288                    11                   11                      4
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 14934_NA                0.0015340846                   0.002438883             10.229588       1.759320e-01
# 20568_NA                0.0411899509                   0.029736338             28.204935       2.018603e-04
# 21991_NA                0.0000000000                   0.001241038              6.413222       4.924104e-01
# 246086_NA               0.0005005612                   0.002324552             15.006678       3.591389e-02
# 66166_NA                0.0161045864                   0.013909035             32.079724       3.925669e-05
# 97848_NA                0.0212018941                   0.015518967             44.849505       1.462740e-07
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 14934_NA                        34.661416                           1.530489               2.2691329
# 20568_NA                        86.939768                           7.656832               1.0058952
# 21991_NA                         0.000000                           1.000000               4.2802655
# 246086_NA                        8.507213                           1.092982               0.7934615
# 66166_NA                        72.469902                           3.632388               1.8441480
# 97848_NA                        82.849838                           5.830849               1.3766620
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained
# 14934_NA                      0.5184600                             0
# 20568_NA                      0.7998255                             0
# 21991_NA                      0.2327476                           100
# 246086_NA                     0.8510303                             0
# 66166_NA                      0.6053753                             0
# 97848_NA                      0.7110143                             0


write.csv(metaOutput, "metaOutput_wCovariates_ADType_Dissection_DepressionModel.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wCovariates_ADType_Dissection_DepressionModel.csv")

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
  
  
  #For DepressionModel:
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,12], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutputFDR, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[36]<-"DepressionModel_FDR"
  
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

# [1] "metaOutputFDR:"
# num [1:16494, 1:36] -0.04776 -0.05277 -0.05996 -0.03335 0.00299 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16494] "23825_114087" "18585_191569" "66514_246307" "20480_65041" ...
# ..$ : chr [1:36] "Log2FC_AD_vs_Ctrl" "Log2FC_ADType_NonTradvsTrad" "Log2FC_Dissection_DG_vs_HPC" "Log2FC_DepressionModel_Yes_vs_No" ...
# NULL
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 181
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_AD_vs_Ctrl Log2FC_ADType_NonTradvsTrad
# 16653_24525               24525               16653        0.09526207                  0.06230221
# 140481_308757            308757              140481       -0.08477753                 -0.04203072
# 70021_290558             290558               70021       -0.24247023                  0.05801205
# 73467_316406             316406               73467        0.10177663                  0.04621843
# 69961_298002             298002               69961       -0.07994736                 -0.04841664
# 71914_305633             305633               71914       -0.35262387                 -0.14353146
# Log2FC_Dissection_DG_vs_HPC Log2FC_DepressionModel_Yes_vs_No      SE_AD  SE_ADType SE_Dissection
# 16653_24525                     0.2792129                      0.017951310 0.01303421 0.02846503    0.02882628
# 140481_308757                  -0.1805559                      0.028663706 0.01357993 0.02411869    0.02440951
# 70021_290558                   -0.4218610                      0.005975652 0.04016353 0.09348114    0.08424272
# 73467_316406                    0.2297901                     -0.037064362 0.01702485 0.03616700    0.03633440
# 69961_298002                   -0.1055425                     -0.022952434 0.01392001 0.03343751    0.02867027
# 71914_305633                   -0.8221705                     -0.183811704 0.06280563 0.12336038    0.10594461
# SE_DepressionModel      pval_AD pval_ADType pval_Dissection pval_DepressionModel    CI_lb_AD
# 16653_24525           0.02590314 2.698989e-13  0.02861664    3.456235e-22            0.4882991  0.06971550
# 140481_308757         0.02686085 4.296576e-10  0.08139266    1.393484e-13            0.2859185 -0.11139371
# 70021_290558          0.08627328 1.569336e-09  0.53487936    5.508841e-07            0.9447793 -0.32118931
# 73467_316406          0.03371422 2.257251e-09  0.20127859    2.543629e-10            0.2716072  0.06840854
# 69961_298002          0.03298004 9.282577e-09  0.14762428    2.320912e-04            0.4864606 -0.10723006
# 71914_305633          0.12656297 1.971019e-08  0.24462123    8.467531e-15            0.1464088 -0.47572065
# CI_lb_ADType CI_lb_Dissection CI_lb_DepressionModel    CI_ub_AD CI_ub_ADType CI_ub_Dissection
# 16653_24525    0.006511765        0.2227144           -0.03281792  0.12080865  0.118092647       0.33571134
# 140481_308757 -0.089302486       -0.2283976           -0.02398260 -0.05816136  0.005241037      -0.13271414
# 70021_290558  -0.125207625       -0.5869737           -0.16311687 -0.16375115  0.241231727      -0.25674835
# 73467_316406  -0.024667580        0.1585760           -0.10314302  0.13514472  0.117104437       0.30100419
# 69961_298002  -0.113952952       -0.1617352           -0.08759213 -0.05266465  0.017119663      -0.04934981
# 71914_305633  -0.385313357       -1.0298181           -0.43187057 -0.22952708  0.098250432      -0.61452284
# CI_ub_DepressionModel Number_Of_Comparisons Number_of_ Contrasts Number_of_Coefficients
# 16653_24525              0.06872054                    22                   22                      4
# 140481_308757            0.08131001                    22                   22                      4
# 70021_290558             0.17506817                    20                   20                      4
# 73467_316406             0.02901430                    20                   20                      4
# 69961_298002             0.04168726                    22                   22                      4
# 71914_305633             0.06424716                    12                   12                      4
# tau2_ResidualHeterogeneity SE_tau2_ResidualHeterogeneity QE_CochransQ_Teststat QEp_CochransQ_pval
# 16653_24525                 0.000000e+00                  0.0006869335              17.24605        0.506254135
# 140481_308757               0.000000e+00                  0.0007307443              15.13402        0.652745370
# 70021_290558                8.924143e-03                  0.0068444867              38.30450        0.001369112
# 73467_316406                0.000000e+00                  0.0014615995               5.15197        0.994944761
# 69961_298002                6.636168e-06                  0.0009042784              22.18315        0.223955260
# 71914_305633                5.822121e-03                  0.0114329210              11.74083        0.163147187
# I2_PercentVar_TrueHeterogeneity H2_Ratio_EffectHetero_overSamplVar QM_ModeratorOmnibusTest
# 16653_24525                         0.0000000                           1.000000               103.01473
# 140481_308757                       0.0000000                           1.000000                54.98991
# 70021_290558                       54.6202097                           2.203624                36.41228
# 73467_316406                        0.0000000                           1.000000                41.33498
# 69961_298002                        0.1818327                           1.001822                14.83141
# 71914_305633                       25.1808291                           1.336556                67.87745
# QMp_Pval_ModeratorOmnibusTest R2_ModeratorVarianceExplained       AD_FDR ADType_FDR Dissection_FDR
# 16653_24525                    3.492854e-22                     100.00000 4.451712e-09  0.5666301   5.700715e-18
# 140481_308757                  6.900293e-12                     100.00000 3.543386e-06  0.7039778   2.873015e-10
# 70021_290558                   6.126432e-08                      82.64334 8.628210e-06  0.9689088   1.081700e-04
# 73467_316406                   5.552285e-09                     100.00000 9.307776e-06  0.8463788   2.097731e-07
# 69961_298002                   1.966538e-03                      99.37859 3.062136e-05  0.8026416   8.413435e-03
# 71914_305633                   1.215249e-14                      95.92982 5.418330e-05  0.8819289   2.327724e-11
# DepressionModel_FDR MouseVsRat_EntrezGene.ID  Mouse_Symbol Mouse_Genetic.Location
# 16653_24525                     1              16653_24525          Kras               Chr6  cM
# 140481_308757                   1            140481_308757        Man2a2               Chr7  cM
# 70021_290558                    1             70021_290558        Nt5dc2              Chr14  cM
# 73467_316406                    1             73467_316406 1700066M21Rik               Chr1  cM
# 69961_298002                    1             69961_298002        Rpp25l               Chr4  cM
# 71914_305633                    1             71914_305633        Antxr2               Chr5  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.                                 Mouse_Name
# 16653_24525                              Chr6:145162425-145195965(-) Kirsten rat sarcoma viral oncogene homolog
# 140481_308757                              Chr7:79998845-80021123(-)                     mannosidase 2, alpha 2
# 70021_290558                              Chr14:30853046-30861081(+)        5'-nucleotidase domain containing 2
# 73467_316406                               Chr1:57416779-57424582(+)                 RIKEN cDNA 1700066M21 gene
# 69961_298002                               Chr4:41712033-41713517(-)         ribonuclease P/MRP 25 subunit-like
# 71914_305633                               Chr5:98032547-98178876(-)                   anthrax toxin receptor 2
#               Rat_Symbol Rat_Genetic.Location Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.
# 16653_24525         Kras             Chr4 q44                                                   NA
# 140481_308757     Man2a2             Chr1 q31                                                   NA
# 70021_290558      Nt5dc2            Chr16 p16                                                   NA
# 73467_316406   C9h2orf69             Chr9 q31                                                   NA
# 69961_298002      Rpp25l             Chr5 q22                                                   NA
# 71914_305633      Antxr2            Chr14 p22                                                   NA
#                                                          Rat_Name
# 16653_24525                           KRAS proto-oncogene, GTPase
# 140481_308757              mannosidase, alpha, class 2A, member 2
# 70021_290558                  5'-nucleotidase domain containing 2
# 73467_316406  similar to human chromosome 2 open reading frame 69
# 69961_298002                  ribonuclease P/MRP subunit p25 like
# 71914_305633                       ANTXR cell adhesion molecule 2

#####################

#Taking a peek

sum(is.na(metaOutputFDR_OrderbyPval$Log2FC_AD_vs_Ctrl))
#[1] 346
#Not terrible, but definitely some loss


sum(metaOutputFDR_OrderbyPval$AD_FDR<0.05, na.rm=TRUE)
#[1] 181

sum(metaOutputFDR_OrderbyPval$ADType_FDR<0.05, na.rm=TRUE)
#[1] 12

sum(metaOutputFDR_OrderbyPval$Dissection_FDR<0.05, na.rm=TRUE)
#[1] 1149

sum(metaOutputFDR_OrderbyPval$DepressionModel_FDR<0.05, na.rm=TRUE)
#[1] 1

metaOutputFDR_OrderbyPval$Mouse_Symbol[which(metaOutputFDR_OrderbyPval$DepressionModel_FDR<0.05)]
#[1] "Zpbp2"

#hmmm.... if depression model genuinely doesn't seem to matter much, I'm inclined to toss it out
#I am worried about overfitting

pdf("Histogram_pvalues_forDepressionModel.pdf", height=5, width=4)
hist(metaOutputFDR_OrderbyPval$pval_DepressionModel)
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


