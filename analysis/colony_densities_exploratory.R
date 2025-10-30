###### Preliminary analysis for bumble bee density and colony persistence paper
### Started October 29, 2025
### J. Melanson

# Set up environment
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

# first, load in packages
source('src/colony_assignment_functions.R')
source('src/colony_density_analysis_src.R')
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(cowplot)
library(ggplot2)
library(raster)
library(sf)
library(tibble)


###########################
### Load in data
###########################
# raw data
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T, fill = TRUE))
sampleEffort2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022sampledata.csv"), sep = ","))
sampleEffort2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"), sep = ","))
samplePoints = as.data.frame(read.table(paste0(bombus_path, "raw_data/allsamplepoints.csv"), sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

# sibships
mix2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mix2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
imp2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
imp2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")

# load landscape raster/shapefiles
landscape1 = raster(paste(bombus_path, "landscape/rasters/FValley_lc_1res.tif", sep = ""))
fv_points = read_sf(paste(bombus_path, "landscape/fvbombus/fvbombus_points.shp", sep = ""))
landcover = read.csv(paste(bombus_path, "landscape/landcover.csv", sep = ""))

# get landscape metrics for each sample_pt
allmetrics = calculateLandscapeMetrics(landcover.raster = landscape1,
                                       site.shapefile = fv_points,
                                       landcover.classification = landcover,
                                       buffer.sizes = c(500))
saveRDS(allmetrics, "data/landscapemetrics/allmetrics.RDS")
fv_summary = allmetrics[[1]]
all_classes = allmetrics[[2]]
all_classes$value_prop = all_classes$value/all_classes$buffer_area


##################################################
# Check for collinearity between landcover types
##################################################
# landcover PCA associated with higher colony density?

# make wide dataframe
all_classes_wide = pivot_wider(all_classes, 
                               id_cols = "sample_pt", 
                               names_from = lc_type, 
                               values_from = value_prop)
all_classes_wide[is.na(all_classes_wide)] = 0


# make correlation matrix
landscape_corr <- cor(
  column_to_rownames(all_classes_wide, var = "sample_pt")
)

#plot
ggcorrplot(landscape_corr, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, ggtheme = ggplot2::theme_minimal())


##################################################
# Plot landscape metrics by site
##################################################
fv_summary$site = substr(fv_summary$sample_pt, start = 1, stop = 1)
ggplot(fv_summary, aes(x = site, y = prop_seminat_500)) +
  geom_boxplot() +
  theme_minimal()
ggplot(fv_summary, aes(x = site, y = landscape_iji_500)) +
  geom_boxplot() +
  theme_minimal()

##################################################
# Visualize colony and worker densities
##################################################

# for mixtus 2022
mix2022 = filter(mix2022, 
                 notes != "male",
                 year == 2022)

# first get number of individuals and colonies per sample pt
m22_workersummary = mix2022 %>%
  group_by(sample_pt, site) %>%
  summarize(num_workers = n(),
            num_colonies = length(unique(sibshipID)))

surveysummary22 = sampleEffort2022 %>%
  group_by(sample_point, site) %>%
  summarize(survey_effort = 5*n())

m22_workersummary = full_join(m22_workersummary, surveysummary22, by = c("sample_pt" = "sample_point", "site" = "site"))
m22_workersummary[is.na(m22_workersummary)] = 0
m22_workersummary$species = "mixtus"
m22_workersummary$year = 2022

# for mixtus 2023
mix2023 = filter(mix2023, 
                 notes != "male",
                 year == 2023)

# first get number of individuals and colonies per sample pt
m23_workersummary = mix2023 %>%
  group_by(sample_pt, site) %>%
  summarize(num_workers = n(),
            num_colonies = length(unique(sibshipID)))

surveysummary23 = sampleEffort2023 %>%
  group_by(sample_point, site) %>%
  summarize(survey_effort = 5*n())

m23_workersummary = full_join(m23_workersummary, surveysummary23, by = c("sample_pt" = "sample_point", "site" = "site"))
m23_workersummary[is.na(m23_workersummary)] = 0
m23_workersummary$species = "mixtus"
m23_workersummary$year = 2023

# for impatiens 2022
imp2022 = filter(imp2022, 
                 notes != "male",
                 year == 2022)

# first get number of individuals and colonies per sample pt
i22_workersummary = imp2022 %>%
  group_by(sample_pt, site) %>%
  summarize(num_workers = n(),
            num_colonies = length(unique(sibshipID)))

i22_workersummary = full_join(i22_workersummary, surveysummary22, by = c("sample_pt" = "sample_point", "site" = "site"))
i22_workersummary[is.na(i22_workersummary)] = 0
i22_workersummary$species = "impatiens"
i22_workersummary$year = 2022


# for impatiens 2023
imp2023 = filter(imp2023, 
                 notes != "male",
                 year == 2023)

# first get number of individuals and colonies per sample pt
i23_workersummary = imp2023 %>%
  group_by(sample_pt, site) %>%
  summarize(num_workers = n(),
            num_colonies = length(unique(sibshipID)))

i23_workersummary = full_join(i23_workersummary, surveysummary23, by = c("sample_pt" = "sample_point", "site" = "site"))
i23_workersummary[is.na(i23_workersummary)] = 0
i23_workersummary$species = "impatiens"
i23_workersummary$year = 2023



# combine all into single df
all_summary = rbind(m22_workersummary,
                    m23_workersummary,
                    i22_workersummary,
                    i23_workersummary)


# visualize num workers vs num colonies relationship
ggplot(all_summary, aes(x = num_workers, y = num_colonies, colour = site)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Number of workers") +
  ylab("Number of colonies") +
  facet_grid(
    rows = vars(species),
    cols = vars(year),
    scales = "fixed",
    space = "fixed") +
  theme_minimal()

##################################################
# Visualize colony temporal span
##################################################

#add julian date to sample effort data frame
sampleEffort2022$date = paste(sampleEffort2022$day, sampleEffort2022$month, sampleEffort2022$year)
sampleEffort2022$date = gsub(" ", "", sampleEffort2022$date, fixed = TRUE)
sampleEffort2022$date <- as.POSIXlt(sampleEffort2022$date, format = "%d%b%y")
sampleEffort2022$julian_date = sampleEffort2022$date$yday

# join sample effort to specimens
mix2022wd = left_join(mix2022, sampleEffort2022, by = c("sample_pt" = "sample_point",
                                                        "year" = "year",
                                                        "site" = "site",
                                                        "round" = "round"))
mix2022wd$sibshipID <- factor(as.character(mix2022wd$sibshipID))


# order sibships by their first date of occurence
mix2022wd = mix2022wd %>%
  group_by(barcode_id) %>%
  mutate(first_date = min(julian_date)) %>%
  ungroup() %>%
  mutate(sibshipID = reorder(factor(barcode_id), first_date))
mix2022wd = mix2022wd %>%
  group_by(sibshipID) %>%
  mutate(first_date = min(julian_date)) %>%
  ungroup() %>%
  mutate(sibshipID = reorder(factor(sibshipID), first_date))


# plot
ggplot(mix2022wd, aes(x = julian_date, y = sibshipID, group = sibshipID)) +
  geom_line(color = "gray50", linewidth = 1.5) +
  geom_point(size = 1.5, color = "lightblue") +
  geom_point(aes(x = julian_date, y = barcode_id)) +
  scale_y_discrete(breaks = NULL, labels = NULL) +
  labs(x = "Julian date", y = "Siblinship") +
  theme_minimal()

