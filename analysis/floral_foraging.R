###### Fit foraging distance model to real data!
### With floral effects and ~multiple timepoints~
### December 5, 2025
### J. Melanson


##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(sf)
library(stringr)

# Prep workspace
# local
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

# remote
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "~/projects/def-ckremen/melanson/"

# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
veg2022 = read.csv(paste0(bombus_path, "/raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "/raw_data/2023vegetationdata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)

prep_stan_floralforaging = function(sibships1,
                                    sibships2,
                                    effort1,
                                    effort2,
                                    veg1,
                                    veg2,
                                    samplepoints){
  
  # Adjust sibships IDs
  sibships2$sibshipID = sibships2$sibshipID + max(sibships1$sibshipID)
  
  # Get colonies x trap x timepoint matrix
  cxl1 = sibships1 %>%
    distinct(site, sibshipID)
  sereduced1 = effort1[effort1$round %in% sibships1$round,colnames(effort1) %in% c("site", "sample_point", "round", "sample_id")]
  CKT1 = left_join(cxl1, sereduced1, by = "site", relationship = "many-to-many")
  CKT1$year = 2022
  
  cxl2 = sibships2 %>%
    distinct(site, sibshipID)
  sereduced2 = effort2[effort2$round %in% sibships2$round,colnames(effort2) %in% c("site", "sample_point", "round", "sample_id")]
  CKT2 = left_join(cxl2, sereduced2, by = "site", relationship = "many-to-many")
  CKT2$year = 2023
  
  CKT = rbind(CKT1, CKT2)
  CKT$round = as.numeric(CKT$round)
  
  # Get site keys
  sitekey = data.frame(site = c("ED", "HR", "NR", "PM", "SD", "W"),
                       siteid = 1:6)
  
  # Get trap keys
  # here we want a different random intercept for each trap x year combo
  # trap "quality" may differ between year so we treat them as unique traps
  CKT$trapyear = paste0(CKT$sample_point, CKT$year)
  trapkey = data.frame(trapyear = sort(unique(CKT$trapyear)),
                       trap_id = 1:length(unique(CKT$trapyear)))
  
  
  # Get floral abundance estimates
  beeflowers = unique(c(specs2022$active_flower, specs2023$active_flower))
  beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]
  
  columnlist = c("sample_id", beeflowers)
  floral_df_wide = veg1[,colnames(veg1) %in% columnlist]
  floral_df_long1 = floral_df_wide %>%
    pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
    mutate(across(-c(sample_id, flower),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(-c(sample_id, flower),
                  ~ 10^(.-1))) %>%
    group_by(sample_id) %>%
    summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
    mutate(across(-c(sample_id),
                  ~ log10(.)))
  
  floral_df_wide = veg2[,colnames(veg2) %in% columnlist]
  floral_df_long2 = floral_df_wide %>%
    pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
    mutate(across(-c(sample_id, flower),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(-c(sample_id, flower),
                  ~ 10^(.-1))) %>%
    group_by(sample_id) %>%
    summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
    mutate(across(-c(sample_id),
                  ~ log10(.)))
  
  floral_df_long = rbind(floral_df_long1, floral_df_long2)
  
  
  # Get trap coordinates
  colnames(samplepoints)[1:2] = c("sample_pt", "coord")
  samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
  samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
  samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
  samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]
  
  # convert to meter based crs
  traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
  traps_sf_m = st_transform(traps_sf, 900913)
  traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
  colnames(traps_m) = c("sample_point", "site", "trap_x", "trap_y")
  traps_m$site = as.factor(traps_m$site)
  
  # Get Julian dates
  #add julian date to sample effort data frame
  effort1$date = paste(effort1$day, effort1$month, effort1$year)
  effort1$date = gsub(" ", "", effort1$date, fixed = TRUE)
  effort1$date <- as.POSIXlt(effort1$date, format = "%d%b%y")
  effort1$julian_date = effort1$date$yday
  
  effort2$date = paste(effort2$day, effort2$month, effort2$year)
  effort2$date = gsub(" ", "", effort2$date, fixed = TRUE)
  effort2$date <- as.POSIXlt(effort2$date, format = "%d%b%y")
  effort2$julian_date = effort2$date$yday
  
  effort = rbind(effort1[,colnames(effort1) %in% c("sample_id", "julian_date")], effort2[,colnames(effort2) %in% c("sample_id", "julian_date")])
  
  
  
  # Get nonzero counts
  counts1 = sibships1 %>%
    filter(notes != "male") %>%
    group_by(round, sample_pt, sibshipID) %>%
    summarize(count=n()) %>%
    ungroup()
  colnames(counts1) = c("round", "sample_point", "sibshipID", "counts")
  
  counts2 = sibships2 %>%
    filter(notes != "male") %>%
    group_by(round, sample_pt, sibshipID) %>%
    summarize(count=n()) %>%
    ungroup()
  colnames(counts2) = c("round", "sample_point", "sibshipID", "counts")
  
  counts = rbind(counts1, counts2)
  
  
  # Join!
  CKT = CKT %>%
    left_join(sitekey) %>%
    left_join(floral_df_long) %>%
    left_join(traps_m) %>%
    left_join(trapkey) %>%
    left_join(counts) %>%
    left_join(effort) %>%
    arrange(round, sibshipID, sample_point)
  
  # Clean up!
  CKT$counts[is.na(CKT$counts)] = 0 # these are true zeroes
  CKT$floral_abundance[CKT$floral_abundance == -Inf] = 0 #questionable -- not true zeroes
  CKT$floral_abundance[is.na(CKT$floral_abundance)] = 0 #questionable -- not true zeroes
  CKT$trap_x = as.numeric(CKT$trap_x)/1000
  CKT$trap_y = as.numeric(CKT$trap_y)/1000
  
  
  # Get stan data!
  data = list(
    C = length(unique(CKT$sibshipID)),
    K = length(unique(CKT$trapyear)),
    O = nrow(CKT),
    trap_pos = cbind(CKT$trap_x, CKT$trap_y),
    colony_id = CKT$sibshipID,
    trap_id = CKT$trap_id,
    fq = CKT$floral_abundance,
    doy = CKT$julian_date,
    yobs = CKT$counts,
    lower_x = (min(CKT$trap_x) - 5),
    upper_x = (max(CKT$trap_x) + 5),
    lower_y = (min(CKT$trap_y) - 5),
    upper_y = (max(CKT$trap_y) + 5)
  )
  
  return(data)
}



### Get stan data for mixtus and impatiens
mixtus_data = prep_stan_floralforaging(mixtus_sibs2022,
                                       mixtus_sibs2023,
                                       effort2022,
                                       effort2023,
                                       veg2022,
                                       veg2023,
                                       samplepoints)

impatiens_data = prep_stan_floralforaging(impatiens_sibs2022,
                                       impatiens_sibs2023,
                                       effort2022,
                                       effort2023,
                                       veg2022,
                                       veg2023,
                                       samplepoints)


# Get task id
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

stanfile = "models/floral_foraging_zinb.stan"

if (task_id == 1){
  stanFit = stan(file = stanfile,
                 data = mixtus_data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 2000,
                 refresh = 10,
                 verbose = TRUE)
  saveRDS(stanFit, "analysis/foraging_modelfits/mixtusfloralfit.rds")
} else{
  stanFit = stan(file = stanfile,
                 data = impatiens_data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 2000,
                 refresh = 10,
                 verbose = TRUE)
  saveRDS(stanFit, "analysis/foraging_modelfits/impatiensfloralfit.rds")
}



# #Plot the posteriors of some colonies
# 
# plot_list = list()
# numplots = 9
# legends = list()
# 
# for (i in 1:numplots){
#   delta_draws = cbind(x = rstan::extract(stanFit, pars = "delta_x")$delta[, i],
#                       y = rstan::extract(stanFit, pars = "delta_y")$delta[, i])
# 
#   p = ggplot(delta_draws, aes(x = x, y = y)) +
#     geom_density_2d_filled(alpha = 0.8) +
# 
#     #plot trap locations / sizes / quality
#     geom_point(data = CKT[CKT$sibshipID ==i,], aes(x = trap_x, y = trap_y, size = counts, colour = "red")) +
#     scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
# 
#     #miscellaneous
#     labs(title = paste("Colony", i),
#          size = "Number of Captures",
#          level = "Colony Posterior") +
#     guides(colour = "none") +
#     #xlim(c(1220, 1225)) +
#     #ylim(c(455,460)) +
#     coord_equal() +
#     theme_bw()
# 
#   # save legend
#   g <- ggplotGrob(p)
#   legend_index <- which(g$layout$name == "guide-box-right")
#   legend <- g$grobs[[legend_index]]
# 
#   # remove legend from plot
#   p <- p + theme(legend.position = "none")
# 
#   #save plot
#   plot_list[[i]] = p
#   legends[[1]] = legend
# }
# 
# fig = grid.arrange(grobs = plot_list, ncol = 3)
# fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
# ggsave(paste0("analysis/colony_posteriors/foragingmodel_", task_id, ".jpg"), fig, height = 3000, width = 4000, units = "px")