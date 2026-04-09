#Code for Gemma Database Extractions
#Sophia Espinoza
#Edited by Eva Geoghegan March 1st 2025
#Adding in reanalysis 
#Adding in GSE230149 from Gemma - parietal cortex

#Create Working Directory
dir.create("./PFC")

file.rename("~/Desktop/BDA_R_Workspace.RData", "PFC/myfile.RData")


#Load code necessary code packages
library(devtools)
library(plyr)
library(tidyverse)
library(gemma.R)
library(dplyr)

#Code for full list of search results:
AllResults <- gemma.R ::get_datasets(query = "cipramil celexa ciprapine lexapro lexam cipralex prozac prozep sarafem selfemra olena faverin luvox seroxat lustral zoloft brintellix elavil vanatrip asendin norpramin nebril irene adapin sinequan tofranil aventyl allegron pamelor vivctil surmontil khedezla pristiq cymbalta drizalma irenka survector wellbutrin zyban merital alival phenotropil carphedon marplan marplon enerzer nardil parnate alnert celeport amira aurorix clobemix depnil manerix eldepryl emsam selgin remeron ludiomil tecipul tolvon edronax qelbree lucelan metatone normarex serzone desyrel viibryd spravato auvelity valdoxan citalopram escitalopram fluoxetine fluvoxamine paroxetine sertraline vortioxetine amitriptyline amoxapine desipramine doxepin imipramine nortriptyline protriptyline trimipramine desvenlafaxine duloxetine amineptine bupropion nomifensine phenylpiracetam tametraline hydrazine isocarboxazid phenelzine tranylcypromine bifemelane moclobemide selegiline mirtazapine maprotiline setiptiline mianserin reboxetine viloxazine teniloxazine etoperidone lorpiprazole nefazodone trazodone lubazodone vilazodone aptazamine ketamine MDMA ECT TMS tianeptine agomelatine psilocybin brexanolone clomipramine dothiepin esketamine iprindole iproniazid lofepramine roboxetine sulpiride venlafaxine DBS SSRI SNRI TCA NDRI MAOI NRIS SARI SPARI NASSA NRISA unicyclic tricyclic tetracyclic antridepress* \"selective serotonin reuptake inhibitor\" \"electroconvulsive therapy\" \"transmagnetic stimulation\" \"deep brain stimulation\" \"serotonin and norepinephrine reuptake inhibitors\" \"norepinephrine dopamine reuptake inhibitor\" \"monoamine oxidase inhibitors\" \"norepinephrine reuptake inhibitor\" \"noradrenergic and specific serotonergic antidepressant\" exercis* \"running wheel\" \"environmental enrichment\" \"cage enrichment\" " ) %>%
  gemma.R:::get_all_pages()

#Examine how many results our search gives us:
str(AllResults)

#Rename full search terms list
MyQueryTerms<-"cipramil celexa ciprapine lexapro lexam cipralex prozac prozep sarafem selfemra olena faverin luvox seroxat lustral zoloft brintellix elavil vanatrip asendin norpramin nebril irene adapin sinequan tofranil aventyl allegron pamelor vivctil surmontil khedezla pristiq cymbalta drizalma irenka survector wellbutrin zyban merital alival phenotropil carphedon marplan marplon enerzer nardil parnate alnert celeport amira aurorix clobemix depnil manerix eldepryl emsam selgin remeron ludiomil tecipul tolvon edronax qelbree lucelan metatone normarex serzone desyrel viibryd spravato auvelity valdoxan citalopram escitalopram fluoxetine fluvoxamine paroxetine sertraline vortioxetine amitriptyline amoxapine desipramine doxepin imipramine nortriptyline protriptyline trimipramine desvenlafaxine duloxetine amineptine bupropion nomifensine phenylpiracetam tametraline hydrazine isocarboxazid phenelzine tranylcypromine bifemelane moclobemide selegiline mirtazapine maprotiline setiptiline mianserin reboxetine viloxazine teniloxazine etoperidone lorpiprazole nefazodone trazodone lubazodone vilazodone aptazamine ketamine MDMA ECT TMS tianeptine agomelatine psilocybin brexanolone clomipramine dothiepin esketamine iprindole iproniazid lofepramine roboxetine sulpiride venlafaxine DBS SSRI SNRI TCA NDRI MAOI NRIS SARI SPARI NASSA NRISA unicyclic tricyclic tetracyclic antridepress* \"selective serotonin reuptake inhibitor\" \"electroconvulsive therapy\" \"transmagnetic stimulation\" \"deep brain stimulation\" \"serotonin and norepinephrine reuptake inhibitors\" \"norepinephrine dopamine reuptake inhibitor\" \"monoamine oxidase inhibitors\" \"norepinephrine reuptake inhibitor\" \"noradrenergic and specific serotonergic antidepressant\" exercis* \"running wheel\" \"environmental enrichment\" \"cage enrichment\""

#Full list of results with no filter:
result_MyQueryTerms_NoFilter<- gemma.R ::get_datasets(query=MyQueryTerms) %>% 
  gemma.R:::get_all_pages() 

#Next, we filter by different characteristics: 
result_MyQueryTerms_NoFilter$taxon.name

result_MyQueryTerms_NoFilter$experiment.troubled

table(result_MyQueryTerms_NoFilter$experiment.rawData)
#All results marked with a -1 means there is no raw data available

#We can filter this down more by searching for rat/mice datasets
result_MyQueryTerms_RatsMice<- gemma.R ::get_datasets(query=MyQueryTerms, taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages()
#This brings our data sets down to 500 

#Examples of grabbing information from the first column:
result_MyQueryTerms_RatsMice$experiment.shortName
result_MyQueryTerms_RatsMice[,1]

#Example of grabbing info from the first row
result_MyQueryTerms_RatsMice[1,]
result_MyQueryTerms_Filtered<-result_MyQueryTerms_RatsMice[result_MyQueryTerms_RatsMice$experiment.troubled==FALSE,]

#Creating a .csv file makes our data easier to read
write.csv(result_MyQueryTerms_Filtered, "result_MyQueryTerms_Filtered.csv")

###################

#Next we want to narrow it down by brain region
MyResults<-result_MyQueryTerms_Filtered

#Create an empty vector to save our results to:
OrganismPartAnnotations_All<-vector(mode="character")

#How many datasets do we have?
length(MyResults$experiment.shortName)
nrow(MyResults)

#Searching for organism parts: 
for(i in c(1:nrow(MyResults))){
  
  ExperimentName<-MyResults$experiment.shortName[i]
  
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  rm(ExperimentName)
  
  #How many annotations do we have?
  length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
    OrganismPartAnnotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
    rm(ExperimentAnnotations)
    
    OrganismPartAnnotations_All<-c(OrganismPartAnnotations_All, OrganismPartAnnotations)
    rm(OrganismPartAnnotations)
    
    #if there isn't organism part annotation, we just tell the loop to move on to the next dataset entry:
    
  }else{
    rm(ExperimentAnnotations)
  }
}
#Results of this search are stored here:
OrganismPartAnnotations_All

#Organizing the data by organism part so it's easier to examine
table(OrganismPartAnnotations_All)

#Save this table:
write.csv(table(OrganismPartAnnotations_All), "Table_OrganismPartAnnotations.csv")

###################

#Next, we want to filter down by specific brain region
result_MyQueryTerms_RatsMice_CerebralCortex <- gemma.R ::get_datasets(query=MyQueryTerms, filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0000956)', taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages()
str(result_MyQueryTerms_RatsMice_CerebralCortex)
#This narrows down our datasets down to 145

#Filtering for high quality data
result_MyQueryTerms_RatsMice_CerebralCortex_Filtered<-result_MyQueryTerms_RatsMice_CerebralCortex[result_MyQueryTerms_RatsMice_CerebralCortex$experiment.troubled==FALSE,]
str(result_MyQueryTerms_RatsMice_CerebralCortex_Filtered)

#We also want to look at filtering out "external" datasets
table(result_MyQueryTerms_RatsMice_CerebralCortex_Filtered$experiment.rawData)
#All datasets marked with a -1 might not be useful
#This leaves us with 112 useful datasets 

#Put into .csv to make a bit more readable
write.csv(result_MyQueryTerms_RatsMice_CerebralCortex_Filtered, "result_MyQueryTerms_RatsMice_CerebralCortex_Filtered.csv")

MyResults<-result_MyQueryTerms_RatsMice_CerebralCortex_Filtered

#Look through organism parts again
OrganismPartAnnotations_All<-vector(mode="character")

for(i in c(1:length(MyResults$experiment.shortName))){
  ExperimentName<-MyResults$experiment.shortName[i]
  
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  rm(ExperimentName)
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
    OrganismPartAnnotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
    rm(ExperimentAnnotations)
    
    OrganismPartAnnotations_All<-c(OrganismPartAnnotations_All, OrganismPartAnnotations)
    rm(OrganismPartAnnotations)
    
  }else{
    rm(ExperimentAnnotations)
  }
}
table(OrganismPartAnnotations_All)

#Save this information to then sort through by hand:
write.csv(OrganismPartAnnotations_All, "OrganismPartAnnotations_AfterFilteringToRegion.csv")

MyResults<-result_MyQueryTerms_RatsMice_CerebralCortex_Filtered

#Let's make some empty vectors that are the same length as the columns in our results
#These empty vectors will be used to store our annotations while we loop through the rows of datasets:
OrganismParts<-vector(mode="character", length=nrow(MyResults))
CellTypes<-vector(mode="character", length=nrow(MyResults))
DevelopmentalStages<-vector(mode="character", length=nrow(MyResults))
Treatments<-vector(mode="character", length=nrow(MyResults))
Diseases<-vector(mode="character", length=nrow(MyResults))
DiseaseModels<-vector(mode="character", length=nrow(MyResults))
Genotypes<-vector(mode="character", length=nrow(MyResults))
Strains<-vector(mode="character", length=nrow(MyResults))
Sex<-vector(mode="character", length=nrow(MyResults))

#I'm going to loop over all of the rows (row number =i) in my results (i.e., dataset metadata)
#And collect all of this annotation information
#And then format it in a way so that it can be added into my simple dataframe of results
#And then outputted and read easily in a spreadsheet program like excel

for(i in c(1:nrow(MyResults))){
  
  #Pulling out the name for the dataset in a row (row number=i):
  ExperimentName<-MyResults$experiment.shortName[i]
  
  #Accessing the annotations for the dataset:
  ExperimentAnnotations<-get_dataset_annotations(dataset=ExperimentName)
  #The number and type of annotations for the datasets is quite variable
  
  rm(ExperimentName)
  
  #Determining whether there is any annotation for organism part:
  
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"])>0){
    
    #If there is organism part annotation, I'm grabbing it:
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="organism part"]
    
    #And then collapsing that vector of annotations into a single string 
    #that can be easily stashed in a single cell in a data.frame (or Excel spreadsheet) 
    #This will eventually become part of the the row for that dataset in the results
    # e.g., "annotation 1; annotation 2; annotation 3"
    OrganismParts[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
    #If there isn't any annotation for organism part, we move on to the next type of annotation:
  }else{}
  
  #Now grabbing the annotation for cell type in a similar manner: 
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="cell type"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="cell type"]
    
    CellTypes[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for developmental stage in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="developmental stage"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="developmental stage"]
    
    DevelopmentalStages[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for treatment in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="treatment"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="treatment"]
    
    Treatments[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for disease in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="disease"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="disease"]
    
    Diseases[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for disease model in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="Disease model"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="Disease model"]
    
    DiseaseModels[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for genotype in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="genotype"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="genotype"]
    
    Genotypes[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for strain in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="strain"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="strain"]
    
    Strains[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  #Now grabbing the annotation for biological sex in a similar manner:
  if(length(ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="biological sex"])>0){
    
    Annotations<-ExperimentAnnotations$term.name[ExperimentAnnotations$class.name=="biological sex"]
    
    Sex[i]<-paste(Annotations, collapse="; ")
    rm(Annotations)
    
  }else{ }
  
  rm(ExperimentAnnotations)
}

#Adding all of those vectors of annotation to my data.frame of results:
MyResults_Annotated<-cbind.data.frame(MyResults, 
                                      OrganismParts,
                                      CellTypes,
                                      DevelopmentalStages,
                                      Treatments,
                                      Diseases,
                                      DiseaseModels,
                                      Genotypes,
                                      Strains,
                                      Sex)

#very pretty, very useful

#Let's add some empty columns for taking inclusion/exclusion notes too

ManipulationUnrelatedToTopic<-vector(mode="character", length=nrow(MyResults))
IncorrectDevelopmentalStage<-vector(mode="character", length=nrow(MyResults))
NotBulkDissection_ParticularCellTypeOrSubRegion<-vector(mode="character", length=nrow(MyResults))
NotFullTranscriptome_ChipSeq_TRAP_miRNA<-vector(mode="character", length=nrow(MyResults))
MetadataIssues_MissingInfo_NoPub_Retracted_Duplicated<-vector(mode="character", length=nrow(MyResults))

Excluded<-vector(mode="character", length=nrow(MyResults))
WhyExcluded<-vector(mode="character", length=nrow(MyResults))

MyResults_Annotated<-cbind.data.frame(MyResults_Annotated, ManipulationUnrelatedToTopic, IncorrectDevelopmentalStage, NotBulkDissection_ParticularCellTypeOrSubRegion, NotFullTranscriptome_ChipSeq_TRAP_miRNA, MetadataIssues_MissingInfo_NoPub_Retracted_Duplicated, Excluded, WhyExcluded)

#And then write out the results so that we can snoop through them in a spreadsheet program like Excel:
write.csv(MyResults_Annotated, "MyResults_Annotated.csv")

################

#We can download the processed expression data for any particular dataset using this code:
Expression<-gemma.R::get_dataset_processed_expression("GSE28644")
Expression<-gemma.R::get_dataset_processed_expression("GSE93041")
Expression<-gemma.R::get_dataset_processed_expression("GSE86392")
Expression<-gemma.R::get_dataset_processed_expression("GSE76694")
Expression<-gemma.R::get_dataset_processed_expression("GSE81672")
Expression<-gemma.R::get_dataset_processed_expression("GSE84183")
Expression<-gemma.R::get_dataset_processed_expression("GSE118670")
Expression<-gemma.R::get_dataset_processed_expression("GSE150264")
Expression<-gemma.R::get_dataset_processed_expression("GSE168172")
Expression<-gemma.R::get_dataset_processed_expression("GSE129359")
Expression<-gemma.R::get_dataset_processed_expression("GSE45229")
Expression<-gemma.R::get_dataset_processed_expression("GSE230149")
str(Expression)

#The first four columns are row metadata: Probe, GeneSymbol, GeneName, NCBIid
#The rest of the columns are gene expression values for each subject

#You can visualize the distribution of the gene expression data for any particular sample using a histogram:
hist(Expression$AI_iTBS_PCTX2_TER)

hist(as.matrix(Expression[,-c(1:4)]))

#We can make this pretty by adding a title and x-axis label:
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression")

#You can also change the color and scaling:
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression", col="red", cex.axis=1.3, cex.lab=1.3)

#You can save the histogram using "export" in the Plots window
#We can also automatically save the histogram by outputting it as a graphics file:
pdf("HistogramForAI_PCTX2.pdf", height=4, width=4)
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression", col="red")
dev.off()
#That will write out into your working directory

#We can also pull out numeric values summarizing the distribution, e.g.:
min(as.matrix(Expression[,-c(1:4)]))
#[1] -5.8601
median(as.matrix(Expression[,-c(1:4)]))
#[1] -2.1651
max(as.matrix(Expression[,-c(1:4)]))
#[1] 12.312
#Or to get more of an overview:
summary(as.matrix(Expression[,-c(1:4)]))

#If we make boxplots for each subject, we see that the data has been normalized:
boxplot(Expression[,-c(1:4)])
#This dataset has been log2 transformed twice!

#If we reverse the log2 transformation by performing 2^X we get a normal-looking distribution:
hist(2^as.matrix(Expression[,-c(1:4)]))

##################

#First Goal:

#Sometimes it is difficult to determine from the abstract and dataset annotation which variables were *actually manipulated* as part of the experiment. 
#This can be especially true for datasets that were part of much larger studies containing multiple experiments
#Gemma's information about the experimental design for the transcriptional profiling experiment can clarify this. 
#This code pulls out that information and places it in a flat format (dataset=row, column=contrast metadata) that can be easily combined with the spreadsheet file where you are already keeping your dataset inclusion/exclusion notes 

ExperimentIDs<-c("GSE28644", "GSE93041", "GSE86392", "GSE76694", "GSE81672", "GSE84183", "GSE118670", "GSE150264", "GSE168172", "GSE129359", "GSE45229", "GSE230149")

#Next, we'll make some empty vectors to store the information that we are going to collect about each dataset:
FactorInfo<-vector(mode="character", length=length(ExperimentIDs))
BaselineFactorInfo<-vector(mode="character", length=length(ExperimentIDs))
TreatmentFactorInfo<-vector(mode="character", length=length(ExperimentIDs))
#And combine them into a nice dataframe:
FactorInfoDF<-cbind.data.frame(ExperimentIDs, FactorInfo, TreatmentFactorInfo, BaselineFactorInfo)

#We'll then loop over each of the datasets:

for(i in c(1:length(ExperimentIDs))){
  #For each dataset, we will use Gemma's API to access the experimental design info:
  Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
  if(nrow(Design)>0){
    #For the dataset, we'll grab the vector of factor categories included in the design
    #And collapse it to a single entry in our data frame.
    FactorInfoDF[i,2]<-paste(Design$factor.category, collapse="; ")
    #For the dataset, we'll next grab the factors included in each result set for the dataset 
    #To do this, we'll first make a vector to store information
    ExperimentalFactorVector<-vector(mode="character", length=1)
    
    #And then loop over each of the result sets for the dataset:
    for(j in c(1:length(Design$result.ID))){
      
      #And grab the factors for the result set and concatenate them to our vector of all factors included in any result set for the dataset:
      ExperimentalFactorVector<-c(ExperimentalFactorVector,Design$experimental.factors[[j]]$summary)
    }
    
    #Then we'll remove the first (empty) entry in the vector and collapse the vector of factors into a single entry for our data frame:
    FactorInfoDF[i,3]<-paste(ExperimentalFactorVector[-1], collapse="; ")
    
    #For the dataset, we'll next grab the baseline/control/reference for each of the factors included in each result set for the dataset 
    #To do this, we'll first make a vector to store information
    BaselineFactorVector<-vector(mode="character", length=1)
    
    #And then loop over each of the result sets for the dataset:
    for(k in c(1:length(Design$result.ID))){
      
      #And grab the baseline/control/reference definitions for the result set and concatenate them to our vector of all baselines included in any result set for the dataset:
      BaselineFactorVector<-c(BaselineFactorVector,Design$baseline.factors[[k]]$summary)
    }
    
    #Then we'll remove the first (empty) entry in the vector and collapse the vector of baselines into a single entry for our data frame:
    FactorInfoDF[i,4]<-paste(BaselineFactorVector[-1], collapse="; ")
    
    #And clean up our workspace before we start the loop again
    rm(Design, ExperimentalFactorVector, BaselineFactorVector)
    
  }else{
    rm(Design)
  }
}

#You can write out this object and use it to help screen datasets
write.csv(FactorInfoDF, "FactorInfoDF.csv")

#################

#Second Goal: 

#For the meta-analysis, we will be extracting the differential expression results from Gemma. The differential expression results for each dataset may include multiple statistical contrasts (e.g., drug1 vs. vehicle, drug2 vs. vehicle). 
#Each of these contrasts are labeled with a result id and contrast id within the Gemma database. 
#We will need to know which of these ids are relevant to our project goals to easily extract their results.
#We will also need to double-check that these statistical contrasts are set up in a manner that makes sense for our experiments:

ResultSets_toScreen<-data.frame(ExperimentID="NA",ResultSetIDs="NA", ContrastIDs="NA", ExperimentIDs="NA", FactorCategory="NA", ExperimentalFactors="NA", BaselineFactors="NA", Subsetted=FALSE, SubsetBy="NA")

str(ResultSets_toScreen)

#We will then loop over each of the datasets:

for(i in c(1:length(ExperimentIDs))){
  
  #For each dataset, we will use Gemma's API to access the experimental design info:
  Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
  
  if(nrow(Design)>0){
    #Next, we'll make some empty vectors to store the experimental factor and baseline factor information for each result id for the dataset:
    ExperimentalFactors<-vector(mode="character", length(Design$result.ID))
    BaselineFactors<-vector(mode="character", length(Design$result.ID))
    
    #We will then loop over each of the result ids for the dataset:
    for(j in c(1:length(Design$result.ID))){
      
      #And grab the vector of experimental factors associated with that result id
      ExperimentalFactorVector<-Design$experimental.factors[[j]]$summary
      #And collapse that info down to a single entry that will fit in our data.frame
      ExperimentalFactors[j]<-paste(ExperimentalFactorVector, collapse="; ")
      
      #And then grab the vector of baseline/control/reference values associated with that result id
      BaselineFactorVector<-Design$baseline.factors[[j]]$summary
      #And collapse that info down to a single entry that will fit in our data.frame
      BaselineFactors[j]<-paste(BaselineFactorVector, collapse="; ")
    }
    
    #Some of the datasets are subsetted for the differential expression analyses
    #We will make an empty vector to store subset information for each result id
    SubsetBy<-vector(mode="character", length(Design$result.ID))
    
    #Then we will determine whether the dataset is subsetted:
    if(Design$isSubset[1]==TRUE){
      
      #If it is subsetted, we will loop over each result id for the dataset
      for (j in c(1:length(Design$result.ID))){
        
        #And grab the vector of subsetting information
        SubsetByVector<-Design$subsetFactor[[j]]$summary
        
        #And then collapse that information down to a single entry that will fit in our dataframe
        SubsetBy[j]<-paste(SubsetByVector, collapse="; ")
      }  
      
      #if the dataset wasn't subsetted for the differential expression analysis:
    }else{
      #We'll just make a vector of NA values to put in the "Subsetted by" column
      SubsetBy<-rep(NA, length((Design$result.ID)))
    }
    
    #Then we combine all of the information for all of the result sets for the dataset into a dataframe
    ResultSets_ForExperiment<-cbind.data.frame(ExperimentID=rep(ExperimentIDs[i],length(Design$result.ID)),ResultSetIDs=Design$result.ID, ContrastIDs=Design$contrast.ID, ExperimentIDs=Design$experiment.ID, FactorCategory=Design$factor.category, ExperimentalFactors, BaselineFactors, Subsetted=Design$isSubset, SubsetBy)
    
    #And add that information as rows to our data frame including the result set information for all datasets:
    ResultSets_toScreen<-rbind.data.frame(ResultSets_toScreen, ResultSets_ForExperiment)
    
    #Then clean up our space before looping to the next dataset:
    rm(ResultSets_ForExperiment, Design, ExperimentalFactors, BaselineFactors, SubsetBy)
    
  }else{
    rm(Design)
  }
  
}

#When we're done, we'll want to remove the initial (empty) row in our data.frame:
ResultSets_toScreen<-ResultSets_toScreen[-1,]

#We can make some empty vectors that we can use to store screening notes:
Include<-vector(mode="character", length=nrow(ResultSets_toScreen))
WrongBaseline<-vector(mode="character", length=nrow(ResultSets_toScreen))
ResultsNotRegionSpecific<-vector(mode="character", length=nrow(ResultSets_toScreen))
ReAnalyze<-vector(mode="character", length=nrow(ResultSets_toScreen))

#And add them as columns to our dataframe:            
ResultSets_toScreen<-cbind.data.frame(ResultSets_toScreen, Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)

#And then write everything out as a .csv file that we can easily mark up in a spreadsheet program:
write.csv(ResultSets_toScreen, "ResultSets_toScreen.csv")

#################

#Goals:
#1. For each of the datasets that have survived our inclusion/exclusion criteria, we will need to identify the relevant result sets and statistical contrasts that are relevant to our research question
#2. We will download the differential expression results for these statistical contrasts
#Consolidate them to one result per gene
#Then pulling out the useful info, calculating sampling variance
#And putting everything in a format that we can easily use in an effect size meta-analysis

list.files()

MyDatasets_Screened<-read.csv("MyDatasets_Screened.csv", stringsAsFactors = FALSE, header=TRUE)

str(MyDatasets_Screened)

#Read in your screened result sets:
setwd("C:/Users/sophi/OneDrive/Documents/BDA_R_Directory")

ResultSet_contrasts<-read.csv("MyDatasets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )
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

#Those names are super unwieldy. I'm going to rename them:
#ComparisonsOfInterest<-c("GSE28644_fluoxetine", "GSE86392_fluoxetine", "GSE76694_ketamine", "GSE81672_ketamine/imipramine", "GSE84183_fluoxetine", "GSE118670_fluoxetine", "GSE150264_imipramine10mg", "GSE168172_sertaline/duloxetine","GSE129359_duloxetine", "GSE45229_quetiapine10mg" )

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

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[3]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[4]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[5]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[6]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[7]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[8]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[9]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[10]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

#######
DE_Results<-differentials[[11]]
str(DE_Results)

FilteringDEResults_GoodAnnotation(DE_Results)

str(DE_Results_GoodAnnotation)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

###############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

source("Function_AligningDEResults.R")

#Aligning the mouse datasets with each other:
ListOfMouseDEResults<-list(DEResults_GSE28644, DEResults_GSE81672, DEResults_GSE118670, DEResults_GSE150264, DEResults_GSE168172, DEResults_GSE129359, DEResults_GSE45229, DEResults_GSE93041)

AligningMouseDatasets(ListOfMouseDEResults)

ListOfRatDEResults<-list(DEResults_GSE76694)

AligningRatDatasets(ListOfRatDEResults)

##################
#Code for aligning the rat and mice results:

MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

#We want to join this ortholog database to our mouse results (Log2FC and SV):
Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)

#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)

#*If there aren't any rat datasets*, we just rename the dataframes so that our downstream code works:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)

MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs
str(MetaAnalysis_SV)

MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")

MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))

##############
#Running the Meta-analysis

#This code caculates the number of NAs (i.e., the number of statistical contrasts lacking real differential expression results) in each row (i.e., for each gene):
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

table(MetaAnalysis_FoldChanges_NAsPerRow)
#I have 11 contrasts and 2 NAs are too many 
NumberOfComparisons=11
CutOffForNAs=7

##########
#Running a basic meta-analysis:
source("Function_RunBasicMetaAnalysis.R")      

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  

str(metaOutput)
head(metaOutput)
tail(metaOutput)

###############
#Accessing some additional gene annotation:

#Reading in a database containing more detailed gene annotation:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat_20240425.csv", header = TRUE, row.names = 1)

colnames(HOM_MouseVsRat)
#Renaming the columns so that we can easily join the annotation to our meta-analysis results:
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)
HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

##################
# Multiple comparison corrections
#This code runs a function that corrects the meta-analysis output to take into account the fact that we are running the statistical calculations thousands of times and therefore have a heightened risk of false discovery (false discovery rate correction) 

source("Function_FalseDiscoveryCorrection.R")

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)

#Make sure that multtest and plyr are loaded or you may come across an error
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")
library(plyr)

#Next we will determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

#Here are the top results as listed by mouse gene symbol:
metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:20)]

#Or as listed by mouse entrez id:
metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:20)]

#Let's plot some of those top results!
hist(metaOutputFDR[,1], breaks=40)

source("Function_MakeForestPlots.R")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="15129", species="Mouse") #Hbb-b1
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="30785", species="Mouse") #Cttnbp2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="241324", species="Mouse") #Crb2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="11966", species="Mouse") #Atp6v1b2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="16369", species="Mouse") #Irs3
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="15275", species="Mouse") #Hk1
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="20965", species="Mouse") #Syn2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="22143", species="Mouse") #Tuba1b
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="230824", species="Mouse") #Grhl3
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="210719", species="Mouse") #Mkx
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="67712", species="Mouse") #Slc25a37
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="13035", species="Mouse") #Ctsg
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="380608", species="Mouse") #Tagap1
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="229214", species="Mouse") #Qrfpr

