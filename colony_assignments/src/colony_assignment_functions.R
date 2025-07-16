##############################################################################
## Custom functions for use in sibship assignment scripts
## Started by J Melanson
## July 14, 2025
##############################################################################

# this function will join Julian dates from sample effort dataframe to the specimen data frame
#add julian date to sample effort data frame
joinJulianDates <- function(specimens, sampleEffort){
  sampleEffort$date = paste(sampleEffort$day, sampleEffort$month, sampleEffort$year)
  sampleEffort$date = gsub(" ", "", sampleEffort$date, fixed = TRUE)
  sampleEffort$date <- as.POSIXlt(sampleEffort$date, format = "%d%b%y")
  justDates = sampleEffort[,colnames(sampleEffort) %in% c("sample_id", "date")]
  joinedSpec = left_join(specimens, justDates, by = "sample_id")
  joinedSpec$julian_date = joinedSpec$date$yday
  return(joinedSpec)
}

# check if something is not in a list/vector
'%!in%' <- function(x,y)!('%in%'(x,y))


# collapse sib pairs into an ordered string
collapse <- function(df) {
  apply(df, 1, function(row) paste(sort(row[c("OffspringID1", "OffspringID2")]), collapse = "-"))
}

# check whether putative families are cliques or not
is_clique <- function(g) {
  v <- vcount(g)
  e <- ecount(g)
  e == v * (v - 1) / 2
}