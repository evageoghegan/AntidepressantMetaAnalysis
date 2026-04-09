#Pulling out validated genes
#Adapted from Dr. Megan Hagenauer
#July 23rd 2025
setwd("/Users/evageoghegan/Brain Data Alchemy Project/GSEA")
list.files()

fGSEA_Results<-read.csv("Non_directional_Sig_Results.csv", header=TRUE, stringsAsFactors = FALSE)

str(fGSEA_Results)
#'data.frame':	10975 obs. of  9 variables:

#Thresholding the fGSEA results by FDR:
fGSEA_Sig<-fGSEA_Results[fGSEA_Results$padj<0.05,]
nrow(fGSEA_Sig)

#Reading in some sort of data frame of DEGs from your meta-analysis:
MetaFDR05_wValidInfo_Exploratory_fGSEA <- read.csv("HPC_whole_DEGS.csv", header = TRUE, stringsAsFactors = FALSE)



str(MetaFDR05_wValidInfo_Exploratory_fGSEA)
#'data.frame':	58 obs. of  11 variables:

colnames(MetaFDR05_wValidInfo_Exploratory_fGSEA)



for (j in 1:length(fGSEA_Sig$leadingEdge)) {
  # Initialize vector to collect matching gene symbols
  MatchingGenes <- c()
  LeadingEdgeVector<-strsplit(fGSEA_Sig$leadingEdge[j], ",")
  print(LeadingEdgeVector[[1]])
  for (i in 1:nrow(MetaFDR05_wValidInfo_Exploratory_fGSEA)) {
    gene <- MetaFDR05_wValidInfo_Exploratory_fGSEA$Mouse_Symbol[i]
    GeneWithAnchors <- paste("^", gene, "$", sep = "")
    print(GeneWithAnchors)
    InGeneSet<- grepl(GeneWithAnchors, LeadingEdgeVector[[1]])
    print(InGeneSet)
    
    
    if(sum(InGeneSet>0)){
      MatchingGenes <- c(MatchingGenes, gene)
      print(MatchingGenes)
    }else{}
  }
  
  fGSEA_Sig$MatchedDEGs[j] <- paste(MatchingGenes, collapse = ";")
}


# View the new column
head(fGSEA_Sig$MatchedDEGs)


# Save the updated results
write.csv(fGSEA_Sig, "nondirectional_Bug_Test_fGSEA_Sig_MatchedDEGs.csv")

###Add in to fix bug: 
i=1 
grepl(paste("^",MetaFDR05_wValidInfo_Exploratory_fGSEA$Mouse_Symbol[i],"$", sep=""), fGSEA_Sig$leadingEdge[j])   

###
#LeadingEdgeVector<-strsplit(fGSEA_Sig$leadingEdge[j])
#print(LeadingEdgeVector)
#InGeneSet<-grepl(paste("^",MetaFDR05_wValidInfo_Exploratory_fGSEA$GeneSymbol[i],"$", sep=""), #LeadingEdgeVector)
#print(InGeneSet
##

LeadingEdgeVector<-strsplit(fGSEA_Sig$leadingEdge[j], ",")
print(LeadingEdgeVector)
InGeneSet<-grepl(paste("^",MetaFDR05_wValidInfo_Exploratory_fGSEA$GeneSymbol[i],"$", sep=""), LeadingEdgeVector)
print(InGeneSet)


