## NaN in GNAR models - investigation 
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
library(spdep) # for neighbourhood construction 
library(expp)
library(rlist)
library(forecast)


# load vectors and GNAR objects 
load(file = "Data/RObjects/GNAR.RData")
load(file = "Data/RObjects/igraph.RData")
load(file = "Data/RObjects/population_weight.RData")
load(file = "Data/RObjects/distance_urbanisation.RData")
load(file = "Data/RObjects/county_index.RData")
load(file = "Data/RObjects/coord_urbanisation.RData")

source("functions_paper.R")

# turn off warnings
options(warn = -1)

# set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Pandemic situations -----------------------------------------------------
load(file = "Data/RObjects/data_subsets_pandemic_situations.RData")

# Queen 
BIC_queen_5 <- fit_and_predict_for_restrictions_all_models(net = covid_net_queen_gnar,
                                                           data_list = datasets_list_coarse,
                                                           forecast_window = 5)

BIC_queen_10 <- fit_and_predict_for_restrictions_all_models(net = covid_net_queen_gnar,
                                                            data_list = datasets_list_coarse,
                                                            forecast_window = 10)

BIC_queen_20 <- fit_and_predict_for_restrictions_all_models(net = covid_net_queen_gnar,
                                                            data_list = datasets_list_coarse,
                                                            forecast_window = 20)

# none for forecasting 5 weeks
BIC_queen_5 %>% pull(BIC) %>% is.nan() %>% which()

# NaN for forecasting 10 weeks
BIC_queen_10 %>% filter(is.nan(BIC))

BIC_queen_20 %>% filter(is.nan(BIC))

# construct GNARfit models: forecasting window 10, GNAR(5, (1, 0, 0, 0, 0))
train_window <- dim(datasets_list_coarse[[1]])[1] - 10

model_1 <- GNARfit(vts = datasets_list_coarse[[1]][1:train_window, ], 
                   net = covid_net_queen_gnar,
                   alphaOrder = 5, 
                   betaOrder = c(1, 0, 0, 0, 0), 
                   globalalpha = TRUE
)

BIC(model_1)


# construct GNARfit models: forecasting window 20, GNAR(1, 0)
train_window <- dim(datasets_list_coarse[[1]])[1] - 20

model_3 <- GNARfit(vts = datasets_list_coarse[[1]][1:train_window, ], 
                   net = covid_net_queen_gnar,
                   alphaOrder = 1, 
                   betaOrder = 0, 
                   globalalpha = TRUE
)

BIC(model_3)


# Economic Hub 
BIC_eco_hub_5 <- fit_and_predict_for_restrictions_all_models(net = covid_net_eco_hubs_gnar, 
                                                             numeric_vertices = TRUE, 
                                                             county_index = county_index_eco_hubs, 
                                                             upper_limit = 4, 
                                                             data_list = datasets_list_coarse)

BIC_eco_hub_10 <- fit_and_predict_for_restrictions_all_models(net = covid_net_eco_hubs_gnar, 
                                                              numeric_vertices = TRUE, 
                                                              county_index = county_index_eco_hubs, 
                                                              upper_limit = 4, 
                                                              data_list = datasets_list_coarse,
                                                              forecast_window = 10)
# none for forecasting 5 weeks
BIC_eco_hub_5 %>% pull(BIC) %>% is.nan() %>% which()

# NaN for forecasting 10 weeks
BIC_eco_hub_10 %>% filter(is.nan(BIC))

model_2 <- GNARfit(vts = datasets_list_coarse[[1]][1:train_window, ], 
                   net = covid_net_eco_hubs_gnar,
                   alphaOrder = 5, 
                   betaOrder = c(2, 0, 0, 0, 0), 
                   globalalpha = TRUE
)

BIC(model_2)
