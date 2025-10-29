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


# function for loading and filtering siblingships

filterSibships = function(filename, 
                          specdata,
                          prob_thresh){
  # get dyads from colony output
  rawdyads = read.table(filename, sep = ",")
  colnames(rawdyads) = c("from", "to", "Probability")
  rawdyads = rawdyads[-1,]
  rawdyads$Probability = as.numeric(rawdyads$Probability)
  
  # filter based on probability threshhold
  dyadsfiltered = rawdyads %>% filter(Probability >= prob_thresh)
  
  # add pair names
  rawdyads$pair = apply(
    rawdyads[, c("from", "to")], 
    1, 
    function(x) paste(sort(x), collapse = "-")
  )
  
  # get edge lists for dyads method
  inferred_edges = dyadsfiltered[,1:2]
  
  # make igraph object
  graph = graph_from_data_frame(inferred_edges, directed = FALSE)
  components = decompose(graph)
  
  # remove cliques to get noncircular families
  noncirculars = components[!sapply(components, is_clique)]
  noncircular_graph = do.call(disjoint_union, noncirculars)
  
  # plot
  if(!is.null(noncircular_graph)){
    plot(noncircular_graph)
  } else (print("No noncircularity."))
  
  
  missing_links = lapply(noncirculars, function(comp) {
    # get all possible pairs of nodes (unordered)
    nodes = V(comp)$name
    node_pairs = t(combn(nodes, 2))
    
    # filter out pairs that already have an edge
    missing <- apply(node_pairs, 1, function(pair) {
      !are_adjacent(comp, pair[1], pair[2])
    })
    
    # make collapsed pair names for each missing link
    missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
      paste(sort(pair), collapse = "-")
    })
    
    df = data.frame(
      missing_pair = missing_names,
      stringsAsFactors = FALSE
    )
    return(df)
  })
  
  missing_df = do.call(rbind, missing_links)
  
  
  # Return to colony outputs!
  # Maintain sib pairs that are highly likely or which are missing but have P > 0.95
  dyadsfiltered_new = rawdyads %>% filter(Probability >= prob_thresh |
                                            (pair %in% missing_df$missing_pair & Probability >= 0.95))
  # Get new edge list
  new_inferred_edges = dyadsfiltered_new[,1:2]
  
  
  # Resolve remaining noncircularity
  # Make new graph
  graph = graph_from_data_frame(new_inferred_edges, directed = FALSE)
  components = decompose(graph)
  
  # For each component, keep only the largest clique
  filtered_components = lapply(components, function(comp) {
    if (is_clique(comp)) {
      return(comp)
    } else {
      cliques = largest_cliques(comp)
      # randomly choose one if there's a tie
      chosen_clique = sample(cliques, 1)[[1]]
      # induce subgraph on that clique
      return(induced_subgraph(comp, chosen_clique))
    }
  })
  
  # Recombine the components into a single graph, make components object
  final_graph = do.call(disjoint_union, filtered_components)
  comp = components(final_graph)
  
  # make dataframe linking each node to its component
  sibship_ids = data.frame(
    barcode_id = V(final_graph)$name,
    sibshipID = comp$membership
  )
  
  # get final dyad edges
  joined = full_join(specdata, sibship_ids, by = "barcode_id") %>%
    mutate(
      sibshipID = if_else(
        is.na(sibshipID),
        # assign new IDs starting from max + 1
        max(sibship_ids$sibshipID, na.rm = TRUE) + 
          cumsum(is.na(sibshipID)),  
        sibshipID
      )
    )
  
  # return final dataset
  return(joined)
}