##### Genotype Simulation Functions #####
# Script Initiated: July 18, 2025
# By: Jenna Melanson
# Goal: functions for simulating:
# (1) genotypes of simulated colonies
# (2) COMIN SOON: single or multiple paternity
# (3) COMING SOON: errors and missing data



simulateGenotypes = function(alleleFreqs,
                            colonyDataDetected
){
  
  # get list of column names (alleles)
  alleles = unique(alleleFreqs$MarkerID)
  df_columns = c()
  for (name in alleles){
    temp = c(paste0(name, "_1"), paste0(name, "_2"))
    df_columns = c(df_columns, temp)
  }
  df_columns = c("individual", "truecolony", df_columns)
  
  #create master dataframe
  sample_df =  setNames(data.frame(matrix(ncol = length(df_columns), nrow = 0)), df_columns)
  
  
  count = 0
  
  for (colony in 1:nrow(colony_data_detected)){
    
    # make a dataframe to hold genotypes for the focal colony
    numsibs = rowSums(yobs_detected)[colony]
    singlesibship_df = as.data.frame(matrix(NA, nrow = numsibs, ncol = length(df_columns)))
    names(singlesibship_df) <- df_columns
    
    # Fill known columns with data
    singlesibship_df$individual = (count+1):(count+numsibs)
    singlesibship_df$truecolony = rep(colony_data_detected$colonyid[colony], numsibs)
    
    #increase count for next sibship df
    count = max(singlesibship_df$individual)
    
    # for each allele, assign some values to the parents and then to each offspring
    allelecounter = 2
    for (allele in unique(impatiens_alellefreq$MarkerID)){
      # assign parental genotypes
      queen = sample(impatiens_alellefreq$AlleleID[impatiens_alellefreq$MarkerID == allele], 
                     size = 2, replace = TRUE, 
                     prob = impatiens_alellefreq$UpdatedFreq[impatiens_alellefreq$MarkerID == allele])
      male = sample(impatiens_alellefreq$AlleleID[impatiens_alellefreq$MarkerID == allele], 
                    size = 1, replace = TRUE, 
                    prob = impatiens_alellefreq$UpdatedFreq[impatiens_alellefreq$MarkerID == allele])
      
      # assign daughter genotypes
      singlesibship_df[,allelecounter + 1] = rep(male, numsibs)
      singlesibship_df[,allelecounter + 2] = sample(queen, size = numsibs, replace = TRUE, prob = c(0.5, 0.5))
      
      allelecounter = allelecounter + 2
    }
    impatienssample_df=rbind(impatienssample_df, singlesibship_df)
  }
  
  
  
}