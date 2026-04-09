#Example usage of Brain.GMT:
#Example Code written for R 4.3.3. using fgsea v.1.2.1 
#Adapted from Toni's
#Changed to run non-directional version

setwd("/Users/evageoghegan/Brain Data Alchemy Project/GSEA")



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.21")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", dependencies = TRUE, ask = FALSE)
BiocManager::install("fgsea", ask = FALSE)



library(clusterProfiler)
library(fgsea)

#Link to the fast gene set enrichment analysis (fGSEA) documentation:
# https://bioconductor.org/packages/release/bioc/html/fgsea.html

#installing the R package fast Gene Set Enrichment Analysis (fGSEA):


#This analysis assumes a differential expression (DE) output file structure similar to that produced by the Limma or EdgeR pipelines 
#Rows=all genes included in the DE analysis, columns=gene annotation and DE statistical output
#At least one of the annotation columns must be official gene symbol
#At least one of the columns of differential statistics must include DE effect size (e.g., Log2 Fold Change)

#Read in the full DE results for a condition from the working directory 
#Replace "DEResults.csv" in the code with your file name
DEResults<-read.csv("metaOutputFDR_annotated.csv", header=TRUE, stringsAsFactors = FALSE)

#Remove rows of DE results that are missing gene symbol annotation or effect size information
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_noNA<-DEResults[is.na(DEResults$Mouse_Symbol)==FALSE & is.na(DEResults$Log2FC_estimate)==FALSE,]

# Create a new column with absolute Log2 Fold Change
DEResults_noNA$abs_Log2FC_estimate <- abs(DEResults_noNA$Log2FC_estimate)

# Average absolute effect sizes per gene symbol
DEResults_absLog2FC_forGSEA <- tapply(X=DEResults_noNA$abs_Log2FC_estimate, INDEX=DEResults_noNA$Mouse_Symbol, FUN=mean)
names(DEResults_absLog2FC_forGSEA) <- names(table(DEResults_noNA$Mouse_Symbol))

# Rank the absolute effect sizes (descending order to prioritize largest changes)
DEResults_absLog2FC_forGSEA_Ranked <- sort(DEResults_absLog2FC_forGSEA, decreasing = TRUE)


#Read in Brain.GMT for your species of interest (this example uses rat)
#If you get a warning about an incomplete line in the .gmt file, just ignore it
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}
library(fgsea)

BrainGMT <- gmtPathways("BrainGMTv2_wGO_MouseOrthologs.gmt.txt")

#Run fast fGSEA on your ranked, averaged effect sizes:
#This code should be compatible with updated fgsea packages - if you have an updated package, this code will run as fgseaSimple()
GSEA_Results <- fgsea(BrainGMT, DEResults_absLog2FC_forGSEA_Ranked, nperm = 10000, minSize = 10, maxSize = 1000)


#Warning message:
#  In readLines(gmt.file) :
#  incomplete final line found on 'BrainGMTv2_wGO_MouseOrthologs.gmt.txt'
#> GSEA_Results <- fgsea(BrainGMT, DEResults_absLog2FC_forGSEA_Ranked, nperm = 10000, minSize = 10, maxSize = 1000)
#Warning messages:
#  1: In fgsea(BrainGMT, DEResults_absLog2FC_forGSEA_Ranked, nperm = 10000,  :
#                You are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To run fgseaMultilevel, you need to remove the nperm #argument in the fgsea function call.
#              2: In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                              All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".
#                                            3: In fgseaSimple(pathways = pathways, stats = stats, minSize = minSize,  :
#                                                              There were 1 pathways for which P-values were not calculated properly due to #unbalanced gene-level statistic values

#Pull out the names for the genes that are driving the enrichment of differential expression in each gene set:
GSEA_Results$leadingEdge<-vapply(GSEA_Results$leadingEdge, paste, collapse= ",", character(1L))

#Write out the results:
write.csv(GSEA_Results, "Non_directional_WholeHPC_BrainGMT_GSEA_Results.csv")

#You can easily view these results in Excel
# Sort by p-value
# padj: false discovery rate (FDR) corrected p-value. This value is normally used to set the threshold for significance (FDR<0.05) 
# ES & NES: Enrichment Score and Normalized Enrichment Score for each gene set. 
# Positive ES & NES values mean that the gene set is enriched with upregulation in response to your variable of interest
# Negative ES & NES values mean that the gene set is enriched with downregulation in response to your variable of interest

# Other aspects of the output can be deciphered by referencing the original GSEA publication: Subramanian et al. 2005
# https://www.pnas.org/doi/10.1073/pnas.0506580102