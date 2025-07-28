##### Genotype Simulation Functions #####
# Script Initiated: July 18, 2025
# By: Jenna Melanson
# Goal: functions for simulating:
# (1) genotypes of simulated colonies
# (2) single or multiple paternity
# (3) COMING SOON: errors and missing data



simulateGenotypes = function(alleleFreqs,
                            colonyDataDetected,
                            observationMatrixDetected,
                            probMultiplePaternity = NULL
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
  for (colony in 1:nrow(colonyDataDetected)){
    
    # make an empty dataframe to hold genotypes for the focal colony
    numsibs = rowSums(observationMatrixDetected)[colony]
    singlesibship_df = as.data.frame(matrix(NA, nrow = numsibs, ncol = length(df_columns)))
    names(singlesibship_df) <- df_columns
    
    # make an empty dataframe for the father genotype (if multiple paternity)
    if (runif(1,0,1) < probMultiplePaternity){
      male_df = as.data.frame(matrix(NA, nrow = 2, ncol = length(alleles)))
      names(male_df) <- alleles
    } else { male_df = NULL}
    
    # Fill known columns with data
    singlesibship_df$individual = (count+1):(count+numsibs)
    singlesibship_df$truecolony = rep(colonyDataDetected$colonyid[colony], numsibs)
    
    #increase count for next sibship df
    count = max(singlesibship_df$individual)
    
    # for each allele, assign some values to the queen
    allelecounter = 2
    
    for (allele in alleles){
      # assign queen genotypes
      queen = sample(alleleFreqs$AlleleID[alleleFreqs$MarkerID == allele], 
                     size = 2, replace = TRUE, 
                     prob = alleleFreqs$UpdatedFreq[alleleFreqs$MarkerID == allele])
      # assign daughter genotypes from mum
      singlesibship_df[,allelecounter + 1] = sample(queen, size = numsibs, replace = TRUE, prob = c(0.5, 0.5))
      
      # assign male genotypes
      # first, determine if the colony will receive multiple fathers
      if (is.null(male_df)){ # simple case--single father
        male = sample(alleleFreqs$AlleleID[alleleFreqs$MarkerID == allele], 
                      size = 1, replace = TRUE, 
                      prob = alleleFreqs$UpdatedFreq[alleleFreqs$MarkerID == allele])
        #assign paternal allele
        singlesibship_df[,allelecounter + 2] = rep(male, numsibs)
        
      } else{ # multiple paternity
        # for simplicity, the maximum number of fathers is 2 -- most colonies have <3 individuals
        #so the impact of additional fathers on observed genotypes would be minimal
        male_df[[allele]] = sample(alleleFreqs$AlleleID[alleleFreqs$MarkerID == allele], 
                           size = 2, replace = TRUE, 
                           prob = alleleFreqs$UpdatedFreq[alleleFreqs$MarkerID == allele])
      }
      
      
      allelecounter = allelecounter + 2
    }
    # assign male allele for multiple paternity case
    if (!is.null(male_df)){
      # pick father for each sib
      rows = sample(c(1,2), size = numsibs, replace = TRUE, prob = c(0.7, 0.3))
      
      for (allele in alleles) {
        col_2 = paste0(allele, "_2")  # corresponding _2 column in sibship df
        # replace each row in the _2 column with the corresponding male allele
          singlesibship_df[[col_2]] = male_df[rows, allele]
        }
      }
    
    sample_df=rbind(sample_df, singlesibship_df)
  }
  sample_df$individual = paste0("Ind", sample_df$individual)
  return(sample_df)
}


# Function to induce realistic error rates and missing data to sibship genotype table
induceErrors = function(genotypeDF,
                        errorRates,
                        missingRates
){
  
  
  
}