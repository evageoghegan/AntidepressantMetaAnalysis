#This code contains a function for making nice forest plots illustrating the effect sizes (Log2FC) for #each of the statistical contrasts and datasets included in our meta-analysis for a particular gene (e.., one of our top findings)
#Megan Hagenauer
#July 25 2024

##########################

#library(metafor)
#Will take away up update all/some/none option

##########################

#Function:

MakeForestPlots <- function(metaOutputFDR_annotated, EntrezIDAsCharacter, species) {
  
  # Conditional for Mouse or Rat input
  if (species == "Mouse") {
    MouseGeneSymbol <- metaOutputFDR_annotated$Mouse_Symbol[which(metaOutputFDR_annotated$Mouse_EntrezGene.ID == EntrezIDAsCharacter)]
    RatGeneSymbol <- metaOutputFDR_annotated$Rat_Symbol[which(metaOutputFDR_annotated$Mouse_EntrezGene.ID == EntrezIDAsCharacter)]
    effect <- as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID == EntrezIDAsCharacter), -c(1:3)])
    var <- as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID == EntrezIDAsCharacter), -c(1:3)]) 
  } else if (species == "Rat") {
    RatGeneSymbol <- metaOutputFDR_annotated$Rat_Symbol[which(metaOutputFDR_annotated$Rat_EntrezGene.ID == EntrezIDAsCharacter)]
    MouseGeneSymbol <- metaOutputFDR_annotated$Mouse_Symbol[which(metaOutputFDR_annotated$Rat_EntrezGene.ID == EntrezIDAsCharacter)]
    effect <- as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID == EntrezIDAsCharacter), -c(1:3)])
    var <- as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID == EntrezIDAsCharacter), -c(1:3)]) 
  } else {
    print("Please use either 'Mouse' or 'Rat' to indicate whether you are using annotation for mouse or rat genes")
    return()
  }
  
  #Define dataset names and corresponding tissue types manually
  tissue_mapping <- list(
    "GSE118670_Fluoxetine" = "Dentate Gyrus",
    "GSE26836_Amitriptyline" = "Ammon's Horn",
    "GSE84183_fluoxetine" = "Dentate Gyrus",
    "GSE43261_Ventral_fluoxetine" = "Dentate Gyrus",
    "GSE43261_Dorsal_fluoxetine" = "Dentate Gyrus",
    "GSE81672_ketamine" = "Ammon's Horn",
    "GSE81672_imipramine" = "Ammon's Horn",
    "GSE73798_ketamine" = "Ammon's Horn",
    "GSE63469_venlafaxine_high" = "Ammon's Horn",
    "GSE63469_venlafaxine_low" = "Ammon's Horn",
    "GSE61301_imipramine_yohimbine" = "Ammon's Horn",
    "GSE56028_imipramine" = "Dentate Gyrus",
    "GSE56028_tianeptine" = "Dentate Gyrus",
    "GSE56028_fluoxetine" = "Dentate Gyrus",
    "GSE56028_agomelatine" = "Dentate Gyrus",
    "GSE27532_desipramine" = "Dentate Gyrus",
    "GSE_230149_TMS_intermittent_theta_burst" = "Ammon's Horn",
    "GSE230148_TMS_modifier_theta_burst" = "Ammon's Horn",
    "GSE230148_TMS_1Hz" = "Ammon's Horn",
    "GSE205325_fluoxetine" = "Ammon's Horn",
    "GSE123027_ECT" = "Ammon's Horn",
    "GSE109445_fluoxetine" = "Ammon's Horn"
  )
  
  # Get dataset names
  dataset_names <- colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)]
  
  # Assign colors based on the tissue type
  tissue_colors <- sapply(dataset_names, function(dataset) {
    tissue_type <- tissue_mapping[[dataset]]
    if (is.null(tissue_type)) {
      return("black")  # Default color if no match is found
    } else if (tissue_type == "Ammon's Horn") {
      return("blue")
    } else if (tissue_type == "Dentate Gyrus") {
      return("red")
    } else {
      return("black")  # Default color for unexpected cases
    }
  })
  
  # Check tissue_colors vector
  print("Assigned Colors:")
  print(tissue_colors)
  
  # Create the Forest Plot
  pdf(paste("ForestPlot_Mouse", MouseGeneSymbol, "Rat", RatGeneSymbol, ".pdf", sep = "_"), height = 8, width = 10)
  


  forest.rma(
    rma(effect, var), 
    slab = dataset_names,  
    xlim = c(-2.25, 2.25),
    ylim = c(0.5, length(dataset_names) + 25),
    colout = tissue_colors, 
    cex = 0.6,  # Reduce text size
    rows = seq(length(dataset_names), 1, by = -1),
    pch = 19,  # Use a smaller, filled plotting symbol
    shade = TRUE
  )
  

  
  # This code labels the forest plot with the mouse and rat gene symbols:
  title(main = MouseGeneSymbol, cex.main = 1.5, line = 1)
  
  
  #This closes the connection to the .pdf file, finishing the plot
  dev.off()
  
}

#Note: the way that this function is currently written, I suspect it might potentially throw up an error message if there is more than one set of log2FC associated with an Entrez ID. 
#This would happen if an Entrez ID in one species (e.g., mouse) mapped to more than one Entrez ID in the other species (e.g., rat) at the point that the results from the two species was joined. 
#It should be solvable by just using the EntrezID for the other species (e.g., rat) as the input to the function


#######################

#Example Usage:

#Note - this function currently uses Mouse Entrez ID (NCBI ID) as it's input
#It needs to be formatted as a character (not an integer) to work
#I will need to make a version of this later that accepts rat Entrez ID

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="13170", species="Mouse") #Dbp

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="74772", species="Mouse") #Atp13a2

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="282580", species="Rat") #Stab2




