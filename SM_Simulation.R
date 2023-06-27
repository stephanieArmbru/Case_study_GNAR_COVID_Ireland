# SIMULATION OF DATA AND VERIFICATION OF GNAR MODELS

rm(list=ls())

# set seed to guarantee reproducibility 
set.seed(1234)

# Load libraries and objects ----------------------------------------------
library(readr)
library(igraph)
library(GNAR)
library(MASS) # for box cox 
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 
library(spdep) # for neighborhood construction 
library(expp)


# load vectors, GNAR objects and data subsets 
load(file = "data/RObjects/GNAR.RData")
load(file = "data/RObjects/county_index.RData")

load(file = "data/RObjects/population_weight.RData")
load(file = "data/RObjects/distance_urbanisation.RData")
load(file = "data/RObjects/coord_urbanisation.RData")

load(file = "data/RObjects/data_subsets_regulations.RData")

# load functions 
source("functions_paper.R")

# turn off warnings
options(warn = -1)

# set ggplot theme 
theme_set(theme_bw(base_size = 16))


# Load data ---------------------------------------------------------------
# load pre-processed data
COVID_weekly_data <- read_csv(file = "data/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)

# transform into correct format for GNAR
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# transform into data frame with time column
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 


# Fit models --------------------------------------------------------------
# fit best performing GNAR model for each network and return GNAR model 
# Queen 
model_queen <- fit_and_predict(alpha = 5, 
                               beta = c(2, 1, 1, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_queen_gnar, 
                               old = TRUE, 
                               county_index = county_index_queen, 
                               return_model = TRUE, 
                               forecast_window = 0)

model_queen$mod


# Simulation Test ---------------------------------------------------------
# exemplarily for Queen's contiguity network
# simulate 120 weeks of data
counties_v <- covid_cases %>% colnames()

# coefficients for GNAR model 
coef_model <- data.frame(type = c("dmatalpha1", "dmatbeta1.1", "dmatbeta1.2",
                                  "dmatalpha2", "dmatbeta2.1", 
                                  "dmatalpha3", "dmatbeta3.1", 
                                  "dmatalpha4",  "dmatbeta4.1",
                                  "dmatalpha5", "dmatbeta5.1"), 
                         param = c(0.1, 0.16, 0.27,
                                   0.15, 0.13,
                                   0.2, 0.05,
                                   -0.1, 0.19,
                                   0.1, -0.1))



simulation_test_queen <- simulate_time_series(gnar_object = covid_net_queen_gnar, 
                                              model_coef = coef_model, 
                                              alpha_order = 5, 
                                              beta_order =  c(2, 1, 1, 1, 1), 
                                              initial_mean = 10, 
                                              var_noise = 1, 
                                              county_index = county_index_queen, 
                                              timeframe = 1000)
save(simulation_test_queen, file = "Data/RObjects/Simulation.RData")

# re-compute GNA model coefficients 
coef_test_queen <- reconstruct_coefficients(model_coef = coef_model, 
                                            simulation_df = simulation_test_queen, 
                                            alphaOrder = 5, 
                                            betaOrder = c(2, 1, 1, 1, 1), 
                                            gnar_object = covid_net_queen_gnar)

# for latex 
strCaption <- "Coefficients for the \\code{GNAR} model serving the data simulation 
and their re-computed values after fitting the \\code{GNAR} model to the simulated data 
and their 95\\% confidence interval (CI) for the \\textbf{Queen's contiguity} network; 
the variance for the random error is set to 1."
print(xtable(coef_test_queen,
             digits=2,
             caption=strCaption,
             label="tab:coef_testing_queen", 
             align = c("", "l", "|", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(coef_test_queen)),
                        command = c(paste("\\toprule \n",
                                          "Coefficient & real value & re-computed value [95\\% CI] \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


