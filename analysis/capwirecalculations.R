############ Make estimates of unobserved colony number using CAPWIRE
### Dec 5 2025
### J Melanson

setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")


# what about CAPWIRE?
library(capwire)
library(dplyr)

mix2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mix2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
imp2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
imp2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")

siblist = list(mix2022, mix2023, imp2022, imp2023)

# create param grid table
# grid = expand.grid(sites = unique(mix2022$site),
#                    sppyear = c("mix2022", "mix2023", "imp2022", "imp2023"))
# key = data.frame(sppyear = c("mix2022", "mix2023", "imp2022", "imp2023"),
#                  id = c(1,2,3,4),
#                  year = c(2022, 2023, 2022, 2023))
# grid = left_join(grid,key)
# grid$taskid = 1:nrow(grid)
# grid$lower = NA
# grid$ML = NA
# grid$upper = NA
# write.csv(grid, "analysis/capwiregrid.csv", row.names = FALSE)


# read in grid to get parameter set
capwiregrid = read.csv("analysis/capwiregrid.csv")

# get task id
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# calculate!
forcap = filter(siblist[[capwiregrid$id[capwiregrid$taskid == task_id]]], 
                notes != "male",
                year == capwiregrid$year[capwiregrid$taskid == task_id],
                site == capwiregrid$sites[capwiregrid$taskid == task_id]) %>%
  group_by(sibshipID) %>%
  summarize(capture.class = n()) %>%
  group_by(capture.class) %>%
  summarize(No.Ind = n()) %>%
  ungroup()

res.tirm = fitTirm(data = as.data.frame(forcap), max.pop = 5000)
conf.int = bootstrapCapwire(fit=res.tirm,
                            bootstraps=1000, CI=c(0.025, 0.975))

# reread file and write results
capwiregrid = read.csv("analysis/capwiregrid.csv")
capwiregrid$ML[capwiregrid$taskid == task_id] = res.tirm$ml.pop.size
capwiregrid$lower[capwiregrid$taskid == task_id] = conf.int$conf.int[1]
capwiregrid$upper[capwiregrid$taskid == task_id] = conf.int$conf.int[2]
write.csv(capwiregrid, "analysis/capwiregrid.csv", row.names = FALSE)

capwiregrid = read.csv("analysis/capwiregrid.csv")
