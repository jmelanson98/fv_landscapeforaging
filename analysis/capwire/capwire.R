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


# Create nice output table for poster
colnames(capwiregrid)[colnames(capwiregrid) == "sites"] = "site"

allsibs = rbind(mix2022, mix2023, imp2022, imp2023)
ngenotyped = allsibs %>%
  group_by(final_id, year, site) %>%
  summarize(Nind = n())
ncolonies = allsibs %>%
  group_by(final_id, year, site) %>%
  summarize(Ncol = length(unique(sibshipID)))

table = capwiregrid %>%
  mutate(
    final_id = case_when(
      substr(sppyear, 1, 3) == "mix" ~ "B. mixtus",
      substr(sppyear, 1, 3) == "imp" ~ "B. impatiens",
      TRUE ~ NA_character_)) %>%
  left_join(ngenotyped, by = c("final_id", "year", "site")) %>%
  left_join(ncolonies, by = c("final_id", "year", "site"))
table$propML = table$Ncol/table$ML
table$proplower = table$Ncol/table$upper
table$propupper = table$Ncol/table$lower
table$Nest = paste0(table$ML, " (", floor(table$lower), ",", ceiling(table$upper), ")")
table$PropDetected = paste0(round(table$propML, digits = 2), " (", round(table$proplower, digits = 2), ",", round(table$propupper, digits = 2), ")")
table$Nest[is.na(table$ML)] = "-"
table$PropDetected[is.na(table$ML)] = "-"


outputcolumns = c("final_id", "site", "Nind", "Ncol", "Nest", "PropDetected")
output22 = table[table$year == 2022, colnames(table) %in% outputcolumns]
output23 = table[table$year == 2023, colnames(table) %in% outputcolumns]
output = left_join(output22, output23, by = c("site", "final_id"))
output$final_id = c("\\emph{B. mixtus}", "", "", "", "", "", "\\emph{B. impatiens}", "", "", "", "", "")


yearnames = c(" ", " ", "2022", " ", " ", " ", "2023", " ", " ", " ")
columnnames =  c("Site", "Species", "$N_{ind}$", "$N_{col}$", "$N_{est}$", "$P(detected)$", "$N_{ind}$", "$N_{col}$", "$N_{est}$", "$P(detected)$")
year_df = as.data.frame(t(yearnames))
column_df = as.data.frame(t(columnnames))
colnames(year_df)   = paste0("V", 1:10)
colnames(column_df) = paste0("V", 1:10)
colnames(output) = paste0("V", 1:10)

output = rbind(year_df, column_df, output)

# Get latex format tex
library(knitr)
library(kableExtra)

kable(output, 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE,
      align = "c",
      col.names = NULL,
      caption = "Number of detected colonies and estimated total number of colonies for \\emph{B. mixtus} and \\emph{B. imaptiens} in 2022 \\& 2023. $N_{ind}$: The number of individuals genotyped from a site; $N_{col}$: The number of unique colonies detected; $N_{est}$: The estimated total number of colonies, with 95\\% bootstrap confidence interval, based on the Two-Innate Rates Model from \\emph{capwire} \\parencite{pennellCapwirePackageEstimating2013}.") %>%
  kable_styling(latex_options = c("striped", "hold_position"))

tbl_latex = kable(
  output,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = NULL,
  align = "c"
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

# Wrap in tikzfigure manually
cat(
  "\\begin{tikzfigure}[]\n",
  "\\caption{Number of detected colonies and estimated total number of colonies for \\emph{B. mixtus} and \\emph{B. imaptiens} in 2022 \\& 2023. $N_{ind}$: The number of individuals genotyped from a site; $N_{col}$: The number of unique colonies detected; $N_{est}$: The estimated total number of colonies, with 95\\% bootstrap confidence interval, based on the Two-Innate Rates Model from \\emph{capwire} \\parencite{pennellCapwirePackageEstimating2013}.}\n",
  "\\centering\n",
  tbl_latex,
  "\n\\end{tikzfigure}\n"
)
