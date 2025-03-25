#load packages
library(gemma.R)
library(dplyr)
library(plyr)


#insert search terms
MyQueryTerms<- gemma.R ::get_datasets(query = "cipramil celexa ciprapine lexapro lexam cipralex prozac prozep sarafem selfemra olena faverin luvox seroxat lustral zoloft brintellix elavil vanatrip asendin norpramin nebril irene adapin sinequan tofranil aventyl allegron pamelor vivctil surmontil khedezla pristiq cymbalta drizalma irenka survector wellbutrin zyban merital alival phenotropil carphedon marplan marplon enerzer nardil parnate alnert celeport amira aurorix clobemix depnil manerix eldepryl emsam selgin remeron ludiomil tecipul tolvon edronax qelbree lucelan metatone normarex serzone desyrel viibryd spravato auvelity valdoxan citalopram escitalopram fluoxetine fluvoxamine paroxetine sertraline vortioxetine amitriptyline amoxapine desipramine doxepin imipramine nortriptyline protriptyline trimipramine desvenlafaxine duloxetine amineptine bupropion nomifensine phenylpiracetam tametraline hydrazine isocarboxazid phenelzine tranylcypromine bifemelane moclobemide selegiline mirtazapine maprotiline setiptiline mianserin reboxetine viloxazine teniloxazine etoperidone lorpiprazole nefazodone trazodone lubazodone vilazodone aptazamine ketamine MDMA ECT TMS tianeptine agomelatine psilocybin brexanolone clomipramine dothiepin esketamine iprindole iproniazid lofepramine roboxetine sulpiride venlafaxine DBS SSRI SNRI TCA NDRI MAOI NRIS SARI SPARI NASSA NRISA unicyclic tricyclic tetracyclic antridepress* \"selective serotonin reuptake inhibitor\" \"electroconvulsive therapy\" \"transmagnetic stimulation\" \"deep brain stimulation\" \"serotonin and norepinephrine reuptake inhibitors\" \"norepinephrine dopamine reuptake inhibitor\" \"monoamine oxidase inhibitors\" \"norepinephrine reuptake inhibitor\" \"noradrenergic and specific serotonergic antidepressant\"") %>%
  gemma.R:::get_all_pages()

MyQueryTerms <- "cipramil celexa ciprapine lexapro lexam cipralex prozac prozep sarafem selfemra olena faverin luvox seroxat lustral zoloft brintellix elavil vanatrip asendin norpramin nebril irene adapin sinequan tofranil aventyl allegron pamelor vivctil surmontil khedezla pristiq cymbalta drizalma irenka survector wellbutrin zyban merital alival phenotropil carphedon marplan marplon enerzer nardil parnate alnert celeport amira aurorix clobemix depnil manerix eldepryl emsam selgin remeron ludiomil tecipul tolvon edronax qelbree lucelan metatone normarex serzone desyrel viibryd spravato auvelity valdoxan citalopram escitalopram fluoxetine fluvoxamine paroxetine sertraline vortioxetine amitriptyline amoxapine desipramine doxepin imipramine nortriptyline protriptyline trimipramine desvenlafaxine duloxetine amineptine bupropion nomifensine phenylpiracetam tametraline hydrazine isocarboxazid phenelzine tranylcypromine bifemelane moclobemide selegiline mirtazapine maprotiline setiptiline mianserin reboxetine viloxazine teniloxazine etoperidone lorpiprazole nefazodone trazodone lubazodone vilazodone aptazamine ketamine MDMA ECT TMS tianeptine agomelatine psilocybin brexanolone clomipramine dothiepin esketamine iprindole iproniazid lofepramine roboxetine sulpiride venlafaxine DBS SSRI SNRI TCA NDRI MAOI NRIS SARI SPARI NASSA NRISA unicyclic tricyclic tetracyclic antridepress* \"selective serotonin reuptake inhibitor\" \"electroconvulsive therapy\" \"transmagnetic stimulation\" \"deep brain stimulation\" \"serotonin and norepinephrine reuptake inhibitors\" \"norepinephrine dopamine reuptake inhibitor\" \"monoamine oxidase inhibitors\" \"norepinephrine reuptake inhibitor\" \"noradrenergic and specific serotonergic antidepressant\""


#Conduct the search with no filter on
result_MyQueryTerms_NoFilter<- gemma.R ::get_datasets(query=MyQueryTerms) %>% 
  gemma.R:::get_all_pages() 

result_MyQueryTerms_NoFilter

str(result_MyQueryTerms_NoFilter)

#to narrow down to a certain brain region and taxa -> hippocampus here
result_MyQueryTerms_RatsMice_Hippocampus <- gemma.R ::get_datasets(query=MyQueryTerms, filter = 'allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/UBERON_0002421)', taxa = c("mouse", "rat")) %>% 
  gemma.R:::get_all_pages() 

str(result_MyQueryTerms_RatsMice_Hippocampus)

#filtered by raw data
result_MyQueryTerms_RatsMice_Hippocampus_Filtered<-result_MyQueryTerms_RatsMice_Hippocampus[result_MyQueryTerms_RatsMice_Hippocampus$experiment.troubled==FALSE,]

str(result_MyQueryTerms_RatsMice_Hippocampus_Filtered)
table(result_MyQueryTerms_RatsMice_Hippocampus_Filtered$experiment.rawData)
write.csv(result_MyQueryTerms_RatsMice_Hippocampus_Filtered, "result_MyQueryTerms_RatsMice_Hippocampus_Filtered.csv")
MyResults<-result_MyQueryTerms_RatsMice_Hippocampus_Filtered


#annotations

MyResults<-result_MyQueryTerms_RatsMice_Hippocampus_Filtered

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