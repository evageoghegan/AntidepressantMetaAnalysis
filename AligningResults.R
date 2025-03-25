#Example pipeline for aligning our results across datasets:
#This has been updated to include the three datasets from renanalysis.R
############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

setwd("/Users/evageoghegan/Brain Data Alchemy Project")

############
library(plyr)
library(dplyr)
#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults<-list(DEResults_GSE123027, DEResults_GSE27532, DEResults_GSE63469, DEResults_GSE73798, DEResults_GSE81672, DEResults_GSE84183, DEResults_GSE26836, DEResults_GSE118670, DEResults_GSE43261_Dorsal, DEResults_GSE43261_Ventral)

AligningMouseDatasets(ListOfMouseDEResults)



#Aligning the rat datasets with each other:

#Example Usage:
#We only have one rat dataset, but we still run it through this code to get it in the right format for our next steps:

ListOfRatDEResults<-list(DEResults_GSE109445, DEResults_GSE205325, DEResults_GSE230149, DEResults_GSE230148, DEResults_GSE56028, DEResults_GSE61301)

AligningRatDatasets(ListOfRatDEResults)


############

#Code for aligning the rat and mice results:

#This code isn't nicely functionalized yet
#It also assumes that there are mouse datasets
#It will break if there are only rat datasets - this needs to be fixed.

################

#First: What are gene orthologs?

# Homology refers to biological features including genes and their products that are descended from a feature present in a common ancestor.

# Homologous genes become separated in evolution in two different ways: separation of two populations with the ancestral gene into two species or gene duplication of the ancestral gene within a lineage:

### Genes separated by speciation are called orthologs.
### Genes separated by gene duplication events are called paralogs.

#This definition came from NCBI (https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html)


#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"

MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	28101 obs. of  7 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	28101 obs. of  7 variables:



#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")

MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")


#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:
colnames(MetaAnalysis_FoldChanges)

# [1] "Rat_EntrezGene.ID"                 "Mouse_EntrezGene.ID"              
# [3] "MouseVsRat_EntrezGene.ID"          "GSE126678_LPS_Acute"              
# [5] "GSE126678_LPS_SubchronicPlusAcute" "GSE126678_LPS_Subchronic"         
# [7] "GSE181285_LPS_Acute"               "GSE205325_LPS_Chronic" 

#Example scatter plot comparing two datasets:
plot(MetaAnalysis_FoldChanges$GSE84183_fluoxetine~MetaAnalysis_FoldChanges$GSE81672_ketamine)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

# Install necessary packages if not already installed
install.packages("ComplexHeatmap")
install.packages("circlize")

# Load the required libraries
library(ComplexHeatmap)
library(circlize)

# Generate the correlation matrix
cor_matrix <- cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use = "pairwise.complete.obs", method = "spearman")

# Define the color gradient using colorRamp2()
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Generate the heatmap with the color gradient
Heatmap(cor_matrix,
        name = "Correlation",  # Name for the color key
        col = col_fun,  # Use the defined color gradient
        cluster_rows = TRUE,  # Cluster rows
        cluster_columns = TRUE,  # Cluster columns
        row_names_gp = gpar(fontsize = 8),  # Adjust row label size
        column_names_gp = gpar(fontsize = 8, rot = 45),  # Adjust column label size and rotate
        heatmap_legend_param = list(title = "Spearman Correlation", 
                                    legend_direction = "horizontal", 
                                    legend_position = "bottom"))  # Adjust legend position



