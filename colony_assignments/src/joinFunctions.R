joinJulianDates <- function(specimens, sampleEffort){
  # this function will join Julian dates from sample effort dataframe to the specimen data frame
  #add julian date to sample effort data frame
  sampleEffort$date = paste(sampleEffort$day, sampleEffort$month, sampleEffort$year)
  sampleEffort$date = gsub(" ", "", sampleEffort$date, fixed = TRUE)
  sampleEffort$date <- as.POSIXlt(sampleEffort$date, format = "%d%b%y")
  justDates = sampleEffort[,colnames(sampleEffort) %in% c("sample_id", "date")]
  joinedSpec = left_join(specimens, justDates, by = "sample_id")
  joinedSpec$julian_date = joinedSpec$date$yday
  return(joinedSpec)
}

