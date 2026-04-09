DownloadingDEResults<-function(ResultSet_contrasts){
  
  #Some ResultSets have more than one statistical contrast, so they are present more than once in our data frame, e.g.:
  #ResultSet_contrasts$ResultSetIDs
  #[1] 553805 553805 553805 570552 556647
  
  #To pull down the statistical results, we'll only want the unique result set ids:
  UniqueResultSetIDs<-unique(ResultSet_contrasts$ResultSetIDs)
  
  print("These are the result sets that you identified as being of interest")
  print(UniqueResultSetIDs)
  #553805 570552 556647
  
  differentials <- UniqueResultSetIDs %>% lapply(function(x){
    #   # take the first and only element of the output. the function returns a list 
    #   # because single experiments may have multiple resultSets. Here we use the 
    #   # resultSet argument to directly access the results we need
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  str(differentials)
  #That code worked. Excellent!
  
  # # some datasets might not have all the advertised differential expression results
  # # calculated due to a variety of factors. here we remove the empty differentials
  missing_contrasts <- differentials %>% sapply(nrow) %>% {.==0}
  #[1] FALSE FALSE FALSE
  differentials <<- differentials[!missing_contrasts]
  UniqueResultSetIDs<<-UniqueResultSetIDs[!missing_contrasts]
  
  print("These are the result sets that had differential expression results:")
  print(UniqueResultSetIDs)
  
  print("Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results")
  
  #Within any particular Result Set, there are likely to be some contrasts that we want and others that we don't want
  #For example, a result set might contain a variety of stress interventions
  #And maybe we only want the acute stress contrast results
  
  #We already identified which statistical contrasts we wanted during our screening:
  #This is the object with the specific contrast ids that we want:
  #ResultSet_contrasts$ContrastIDs
  
  #Which will be these columns within the listed dataframes of differential expression results:
  
  print("These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes")
  Contrasts_Log2FC<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_log2fc", sep="")
  
  print(Contrasts_Log2FC)
  #[1] "contrast_151618_log2fc" "contrast_151617_log2fc" "contrast_151619_log2fc"
  #[4] "contrast_186753_log2fc" "contrast_204289_log2fc"
  
  print("these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances")
  
  Contrasts_Tstat<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_tstat", sep="")
  
  print(Contrasts_Tstat)
  # [1] "contrast_151618_tstat" "contrast_151617_tstat" "contrast_151619_tstat"
  # [4] "contrast_186753_tstat" "contrast_204289_tstat"  
  
}
#This function uses unique because sometimes the same results set ID comes up 
#when you hit print(UniqueResultsSetID) you should be able to see the IDs that you are looking for 
#Line 14 for differentials means lapply (run function over again) for get differential results 
#Then gives you a list with dataframes 
#Within Gemma sometimes you will have a results set with no results - issue with the data 
#Then use the code to take away any results sets that have no information in them
#<< moves the object to the next environment above that
#we are interested in log2fc changes - differences in the expression of genes 
#-1 log2fc means halves
#1 log2fc menas doubled 
#then grab t stats as well
#Because then you can calculate the standard error for the analysis
#line 47 - paste is a way to create strings 