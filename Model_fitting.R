### CHANGE IN GNAR MODELS FOR DATA SUBSETS ACCORDING TO COVID-19 
### REGULATIONS 

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

all_counties <- colnames(datasets_list_coarse[[1]])

# ARIMA benchmark ---------------------------------------------------------
arima_results_list <- lapply(datasets_list_coarse, FUN = function(i) {
  results_arima <- list()
  start_date_year <- i %>% 
    rownames() %>% 
    as.Date() %>% 
    min() %>% 
    substr(start = 1, stop = 4) %>% 
    as.numeric()
  
  start_date_month <- i %>% 
    rownames() %>% 
    as.Date() %>% 
    min() %>% 
    substr(start = 6, stop = 7) %>% 
    as.numeric()
  
  
  
  for (county in i %>% colnames()) {
    
    covid_cases_county <- i %>% 
      as.data.frame() %>% 
      dplyr::select(county %>% all_of()) %>% 
      ts(frequency = 52, 
         start = c(start_date_year, start_date_month))
    
    arima_model <- auto.arima(y = covid_cases_county, 
                              d = 0)
    
    results_arima[[county]] <- list("model" = arima_model, 
                                    data.frame("BIC" = arima_model %>% BIC(), 
                                               "AIC" = arima_model %>% AIC()), 
                                    "arma" = paste0(arima_model$arma, 
                                                    collapse = "-"))
  }
  
  results_arima %>% return()
})


mean_arima_BIC <- lapply(arima_results_list, FUN = function(j) {
  lapply(j, FUN = function(i) {
    i[[2]]
  }) %>% 
    list.rbind() %>% 
    summarise(mean_BIC = mean(BIC), 
              mean_AIC = mean(AIC))
}) %>% list.rbind()

# p - q - d
lapply(arima_results_list, FUN = function(j) {
  lapply(j, FUN = function(i) {
    i[[3]]
  }) %>% list.rbind() %>% table()
})

# predict and compute MASE for ARIMA models 
mase_arima_restrictive <- fit_and_predict_arima(forecast_window = 5, 
                                                results = arima_results_list[[1]], 
                                                data = datasets_list_coarse[[1]] %>% 
                                                  as.data.frame(), 
                                                counties = all_counties
                                                )
mase_arima_free <- fit_and_predict_arima(forecast_window = 5, 
                                         results = arima_results_list[[2]], 
                                         data = datasets_list_coarse[[2]] %>% 
                                           as.data.frame(), 
                                         counties = all_counties
                                         )

# Best performing GNAR models ---------------------------------------------
# compute GNAR models for each data subsets and select the best performing one 
# based on the BIC 

# Queen 
best_for_subset_queen <- fit_and_predict_for_restrictions(net = covid_net_queen_gnar,
                                                          data_list = datasets_list_coarse)
best_for_subset_queen$network <- "Queen"

# Economic hub
best_for_subset_eco_hub <- fit_and_predict_for_restrictions(net = covid_net_eco_hubs_gnar, 
                                                            numeric_vertices = TRUE, 
                                                            county_index = county_index_eco_hubs, 
                                                            upper_limit = 4, 
                                                            data_list = datasets_list_coarse)
best_for_subset_eco_hub$network <- "Eco. hub"

# Railway-based 
best_for_subset_train <- fit_and_predict_for_restrictions(net = covid_net_train_gnar, 
                                                          numeric_vertices = TRUE, 
                                                          county_index = county_index_train, 
                                                          data_list = datasets_list_coarse)
best_for_subset_train$network <- "Train"

# Delaunay triangulation 
best_for_subset_delaunay <- fit_and_predict_for_restrictions(net = covid_net_delaunay_gnar, 
                                                             data_list = datasets_list_coarse)
best_for_subset_delaunay$network <- "Delaunay"

# Gabriel 
best_for_subset_gabriel <- fit_and_predict_for_restrictions(net = covid_net_gabriel_gnar, 
                                                            data_list = datasets_list_coarse)
best_for_subset_gabriel$network <- "Gabriel"

# Relative neighbourhood
best_for_subset_relative <- fit_and_predict_for_restrictions(net = covid_net_relative_gnar, 
                                                             data_list = datasets_list_coarse)
best_for_subset_relative$network <- "Relative"

# SOI 
best_for_subset_soi <- fit_and_predict_for_restrictions(net = covid_net_soi_gnar, 
                                                        data_list = datasets_list_coarse)
best_for_subset_soi$network <- "SOI"

# Complete
best_for_subset_complete <- fit_and_predict_for_restrictions(net = complete_net_gnar, 
                                                             upper_limit = 1, 
                                                             data_list = datasets_list_coarse)
best_for_subset_complete$network <- "Complete"

# compare best-performing models across networks
best_for_subset <- rbind(best_for_subset_queen,
                         best_for_subset_eco_hub, 
                         best_for_subset_train,
                         best_for_subset_delaunay, 
                         best_for_subset_gabriel, 
                         best_for_subset_relative, 
                         best_for_subset_soi, 
                         best_for_subset_complete)

# for latex 
strCaption <- "Overview over the best performing model for every COVID-19 data 
subset for every COVID-19 network, excluding the KNN and DNN network"
print(xtable(best_for_subset[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_subsets", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_for_subset[, c(4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Network & data subset & best model &
                                          BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# construct KNN networks with different neighbourhood size k, fit GNAR models 
# and select best performing model via the BIC
knn_best <- list()

for (k in seq(1, 26, by = 2)) {
  # create nb list
  nb_knn <- knearneigh(x = coord_urbanisation, 
                       k = k, 
                       longlat = TRUE) %>% 
    knn2nb(row.names = coord_urbanisation %>% row.names())
  
  # create igraph object
  covid_net_knn_igraph <- neighborsDataFrame(nb = nb_knn) %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    igraph::simplify() 
  
  # create GNAR object 
  covid_net_knn <- covid_net_knn_igraph %>% 
    igraphtoGNAR()
  
  # create ordered county index data frame 
  county_index_knn <- data.frame("CountyName" = covid_net_knn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  
  # compute an upper limit for neighbourhood stage 
  upper_limit_knn <- distance_table(covid_net_knn_igraph)$res %>% length()
  
  # fit GNAR models and select the best performing one for each data subset 
  res <- fit_and_predict_for_restrictions(net = covid_net_knn, 
                                          upper_limit = upper_limit_knn, 
                                          data_list = datasets_list_coarse)
  
  res$hyperparam <- k
  
  # save best performing model for every k across all data subsets  
  knn_best[[length(knn_best) + 1]] <- res
  
}

# filter the best performing GNAR model for each data subset across all 
# neighbourhood sizes 
knn_best_df <- do.call(rbind.data.frame, knn_best) %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  arrange(data_subset)
knn_best_df$network <-  "KNN"



# construct DNN networks for different distance thresholds d, fit GNAR models 
# and select the best performing one via the BIC 
dnn_best <- list()

for (d in seq(100, 
              350, 
              by = 25)) {
  # create nb list
  nb_dnn <- dnearneigh(x = coord_urbanisation, 
                       d1 = 0, 
                       d2 = d,
                       row.names = coord_urbanisation %>% rownames(),
                       longlat = TRUE, 
                       use_s2 = TRUE)
  # create igraph object 
  covid_net_dnn_igraph <- neighborsDataFrame(nb = nb_dnn) %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    igraph::simplify() 
  
  # create GNAR object 
  covid_net_dnn <- covid_net_dnn_igraph %>% 
    igraphtoGNAR()
  
  # create ordered county index data frame 
  county_index_dnn <- data.frame("CountyName" = covid_net_dnn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  
  # compute an upper limit for neighbourhood stage 
  upper_limit_dnn <- distance_table(covid_net_dnn_igraph)$res %>% length()
  
  # fit GNAR models and select the best performing one for each data subset
  res <- fit_and_predict_for_restrictions(net = covid_net_dnn, 
                                          upper_limit = upper_limit_dnn, 
                                          data_list = datasets_list_coarse)
  
  res$hyperparam <- d
  
  # save best performing model for every distance threshold d
  dnn_best[[length(dnn_best) + 1]] <- res
  
}

# filter best performing model across distance threshold for each data subset
dnn_best_df <- do.call(rbind.data.frame, dnn_best) %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  arrange(data_subset)
dnn_best_df$network <- "DNN"


best_subset_knn_dnn <- rbind(knn_best_df,
                             dnn_best_df)

# for latex 
strCaption <- "Overview over the best performing model and optimal 
neighbourhood size $k$ / distance threshold $d$ for the KNN and DNN network"
print(xtable(best_subset_knn_dnn[, c(5, 1, 4, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_knn_dnn_subsets", 
             align = c("", "l", "|", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_subset_knn_dnn[, c(5, 4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Network & data subset & k / d [in km] &
                                          \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# create a data frame with the best performing GNAR model and network for 
# each data subset 
best_subset_knn_dnn_final <- best_subset_knn_dnn %>% 
  mutate(network = paste(network, hyperparam, sep = "-")) %>% 
  dplyr::select(-hyperparam)


best_for_subset_large_df <- rbind(best_for_subset, 
                                  best_subset_knn_dnn_final)

best_for_subset_all <-  best_for_subset_large_df %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  arrange(data_subset)

# Construct optimal network (KNN / DNN) -----------------------------------
# DNN d = 200
dnn_200 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 200,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_200_igraph<- neighborsDataFrame(nb = dnn_200) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_200_gnar <- dnn_200_igraph %>% 
  igraphtoGNAR()

# DNN d = 325
dnn_325 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 325,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_325_igraph<- neighborsDataFrame(nb = dnn_325) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_325_gnar <- dnn_325_igraph %>% 
  igraphtoGNAR()

# KNN k = 7
knn_7 <- knearneigh(x = coord_urbanisation, 
                     k = 7, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_7_igraph<- neighborsDataFrame(nb = knn_7) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_7_gnar <- knn_7_igraph %>% 
  igraphtoGNAR()

# KNN k = 21
knn_21 <- knearneigh(x = coord_urbanisation, 
                     k = 21, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_21_igraph<- neighborsDataFrame(nb = knn_21) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_21_gnar <- knn_21_igraph %>% 
  igraphtoGNAR()

# Data formatting ---------------------------------------------------------
# format data subsets as data frame with time column 
data_1 <- datasets_list_coarse[[1]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list_coarse[[1]]) %>% as.Date())

data_2 <- datasets_list_coarse[[2]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list_coarse[[2]])  %>% as.Date())



# MASE for restrictive ----------------------------------------------------
# fit best performing GNAR model for each network for data subset 1
best_for_subset %>% filter(data_subset == 1)
best_subset_knn_dnn %>% filter(data_subset == 1)

# Queen
mod_1_queen <- fit_and_predict(alpha = 5, 
                               beta = c(2, 1, 1, 1, 1), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list_coarse[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_1_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(2, 1, 1, 1, 1), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list_coarse[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_1_train <- fit_and_predict(alpha = 5, 
                               beta = c(3, 1, 1, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list_coarse[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)

# Delaunay
mod_1_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 0, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list_coarse[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_1_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 1, 1, 1, 1), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list_coarse[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_1_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(5, 1, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list_coarse[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_1_soi <- fit_and_predict(alpha = 5, 
                             beta = c(2, 0, 0, 0, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list_coarse[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_1_knn <- fit_and_predict(alpha = 5, 
                             beta = c(4, 2, 2, 2, 0), 
                             net = knn_7_gnar, 
                             vts = datasets_list_coarse[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_1_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 1), 
                             net = dnn_200_gnar, 
                             vts = datasets_list_coarse[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_1_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 0, 0, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list_coarse[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# compute MASE for best performing GNAR model for each network
mase_1_queen <- compute_MASE(model = mod_1_queen, 
                             network_name = "subset_1_queen", 
                             n_ahead = 5, 
                             data_df = data_1, 
                             counties = all_counties)

mase_1_eco_hub <- compute_MASE(model = mod_1_eco_hub, 
                               network_name = "subset_1_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_1, 
                               counties = all_counties)

mase_1_train <- compute_MASE(model = mod_1_train, 
                             network_name = "subset_1_train", 
                             n_ahead = 5, 
                             data_df = data_1, 
                             counties = all_counties)

mase_1_delaunay <- compute_MASE(model = mod_1_delaunay, 
                                network_name = "subset_1_delaunay", 
                                n_ahead = 5, 
                                data_df = data_1, 
                                counties = all_counties)

mase_1_gabriel <- compute_MASE(model = mod_1_gabriel, 
                               network_name = "subset_1_gabriel", 
                               n_ahead = 5, 
                               data_df = data_1, 
                               counties = all_counties)

mase_1_relative <- compute_MASE(model = mod_1_relative, 
                                network_name = "subset_1_relative", 
                                n_ahead = 5, 
                                data_df = data_1, 
                                counties = all_counties)

mase_1_soi <- compute_MASE(model = mod_1_soi, 
                           network_name = "subset_1_soi", 
                           n_ahead = 5, 
                           data_df = data_1, 
                           counties = all_counties)

mase_1_knn <- compute_MASE(model = mod_1_knn, 
                           network_name = "subset_1_knn", 
                           n_ahead = 5, 
                           data_df = data_1, 
                           counties = all_counties)

mase_1_dnn <- compute_MASE(model = mod_1_dnn, 
                           network_name = "subset_1_dnn", 
                           n_ahead = 5, 
                           data_df = data_1, 
                           counties = all_counties)

mase_1_complete <- compute_MASE(model = mod_1_complete, 
                                network_name = "subset_1_complete", 
                                n_ahead = 5, 
                                data_df = data_1, 
                                counties = all_counties)


mase_1_overview <- rbind(mase_1_queen, 
                         mase_1_eco_hub, 
                         mase_1_train, 
                         mase_1_delaunay, 
                         mase_1_gabriel, 
                         mase_1_relative, 
                         mase_1_soi, 
                         mase_1_knn, 
                         mase_1_dnn, 
                         mase_1_complete, 
                         mase_arima_restrictive)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network

m_1_delaunay_I <- plot_mase_I(mase_overview = mase_1_overview, 
                              counties_subset = all_counties[1:9],
                              number_counties = 1)
m_1_delaunay_II <- plot_mase_I(mase_overview = mase_1_overview, 
                               counties_subset = all_counties[10:18], 
                               number_counties = 2)
m_1_delaunay_III <- plot_mase_I(mase_overview = mase_1_overview, 
                                counties_subset = all_counties[19:26], 
                                number_counties = 3)

m_1_knn_I <- plot_mase_II(mase_overview = mase_1_overview, 
                          counties_subset = all_counties[1:9],
                          number_counties = 1)
m_1_knn_II <- plot_mase_II(mase_overview = mase_1_overview, 
                           counties_subset = all_counties[10:18],
                           number_counties = 2)
m_1_knn_III <- plot_mase_II(mase_overview = mase_1_overview, 
                            counties_subset = all_counties[19:26],
                            number_counties = 3)

# Predicted vs. fitted for restricted -------------------------------------
g_1_delaunay_I <- plot_predicted_vs_fitted_I(mase_overview = mase_1_overview, 
                                             counties_subset = all_counties[1:9],
                                             number_counties = 1)
g_1_delaunay_II <- plot_predicted_vs_fitted_I(mase_overview = mase_1_overview, 
                                              counties_subset = all_counties[10:18], 
                                              number_counties = 2)
g_1_delaunay_III <- plot_predicted_vs_fitted_I(mase_overview = mase_1_overview, 
                                               counties_subset = all_counties[19:26], 
                                               number_counties = 3)

g_1_knn_I <- plot_predicted_vs_fitted_II(mase_overview = mase_1_overview, 
                                         counties_subset = all_counties[1:9],
                                         number_counties = 1)
g_1_knn_II <- plot_predicted_vs_fitted_II(mase_overview = mase_1_overview, 
                                         counties_subset = all_counties[10:18],
                                         number_counties = 2)
g_1_knn_III <- plot_predicted_vs_fitted_II(mase_overview = mase_1_overview, 
                                         counties_subset = all_counties[19:26],
                                         number_counties = 3)


# MASE for free -----------------------------------------------------------
# fit best performing GNAR model for each network for data subset 1
best_for_subset %>% filter(data_subset == 2)
best_subset_knn_dnn %>% filter(data_subset == 2)

# Queen
mod_2_queen <- fit_and_predict(alpha = 5, 
                               beta = c(3, 0, 0, 0, 0), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list_coarse[[2]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_2_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(3, 2, 2, 2, 0), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list_coarse[[2]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_2_train <- fit_and_predict(alpha = 5, 
                               beta = c(5, 0, 0, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list_coarse[[2]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_2_delaunay <- fit_and_predict(alpha = 4, 
                                  beta = c(4, 1, 1, 1), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list_coarse[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_2_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 0, 0, 0, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list_coarse[[2]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_2_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(5, 0, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list_coarse[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_2_soi <- fit_and_predict(alpha = 5, 
                             beta = c(4, 1, 0, 0, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list_coarse[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_2_knn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 0), 
                             net = knn_21_gnar, 
                             vts = datasets_list_coarse[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_2_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 0), 
                             net = dnn_325_gnar, 
                             vts = datasets_list_coarse[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_2_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 1, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list_coarse[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_2_queen <- compute_MASE(model = mod_2_queen, 
                             network_name = "subset_2_queen", 
                             n_ahead = 5, 
                             data_df = data_2, 
                             counties = all_counties)

mase_2_eco_hub <- compute_MASE(model = mod_2_eco_hub, 
                               network_name = "subset_2_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_2, 
                               counties = all_counties)

mase_2_train <- compute_MASE(model = mod_2_train, 
                             network_name = "subset_2_train", 
                             n_ahead = 5, 
                             data_df = data_2,
                             counties = all_counties)

mase_2_delaunay <- compute_MASE(model = mod_2_delaunay, 
                                network_name = "subset_2_delaunay", 
                                n_ahead = 5, 
                                data_df = data_2, 
                                counties = all_counties)

mase_2_gabriel <- compute_MASE(model = mod_2_gabriel, 
                               network_name = "subset_2_gabriel", 
                               n_ahead = 5, 
                               data_df = data_2, 
                               counties = all_counties)

mase_2_relative <- compute_MASE(model = mod_2_relative, 
                                network_name = "subset_2_relative", 
                                n_ahead = 5, 
                                data_df = data_2, 
                                counties = all_counties)

mase_2_soi <- compute_MASE(model = mod_2_soi, 
                           network_name = "subset_2_soi", 
                           n_ahead = 5, 
                           data_df = data_2, 
                           counties = all_counties)

mase_2_knn <- compute_MASE(model = mod_2_knn, 
                           network_name = "subset_2_knn", 
                           n_ahead = 5, 
                           data_df = data_2, 
                           counties = all_counties)

mase_2_dnn <- compute_MASE(model = mod_2_dnn, 
                           network_name = "subset_2_dnn", 
                           n_ahead = 5, 
                           data_df = data_2, 
                           counties = all_counties)

mase_2_complete <- compute_MASE(model = mod_2_complete, 
                                network_name = "subset_2_complete", 
                                n_ahead = 5, 
                                data_df = data_2, 
                                counties = all_counties)


mase_2_overview <- rbind(mase_2_queen, 
                         mase_2_eco_hub, 
                         mase_2_train, 
                         mase_2_delaunay, 
                         mase_2_gabriel, 
                         mase_2_relative, 
                         mase_2_soi, 
                         mase_2_knn, 
                         mase_2_dnn, 
                         mase_2_complete, 
                         mase_arima_free)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network


m_2_delaunay_I <- plot_mase_I(mase_overview = mase_2_overview, 
                              mase_name = "free", 
                              counties_subset = all_counties[1:9],
                              number_counties = 1, 
                              types = c("subset_2_delaunay", 
                                        "subset_2_gabriel", 
                                        "subset_2_relative", 
                                        "subset_2_soi", 
                                        "subset_2_train", 
                                        "ARIMA"), 
                              color_types = c("ARIMA" = "grey", 
                                              "subset_2_gabriel" = "#00BFC4", 
                                              "subset_2_relative" = "#00B0F6", 
                                              "subset_2_soi" = "#9590FF", 
                                              "subset_2_delaunay" = "#E76BF3", 
                                              "subset_2_train" = "#FF62BC"))
m_2_delaunay_II <- plot_mase_I(mase_overview = mase_2_overview,
                               mase_name = "free", 
                               counties_subset = all_counties[10:18], 
                               number_counties = 2, 
                               types = c("subset_2_delaunay", 
                                         "subset_2_gabriel", 
                                         "subset_2_relative", 
                                         "subset_2_soi", 
                                         "subset_2_train", 
                                         "ARIMA"), 
                               color_types = c("ARIMA" = "grey", 
                                               "subset_2_gabriel" = "#00BFC4", 
                                               "subset_2_relative" = "#00B0F6", 
                                               "subset_2_soi" = "#9590FF", 
                                               "subset_2_delaunay" = "#E76BF3", 
                                               "subset_2_train" = "#FF62BC"))
m_2_delaunay_III <- plot_mase_I(mase_overview = mase_2_overview, 
                                mase_name = "free", 
                                counties_subset = all_counties[19:26], 
                                number_counties = 3, 
                                types = c("subset_2_delaunay", 
                                          "subset_2_gabriel", 
                                          "subset_2_relative", 
                                          "subset_2_soi", 
                                          "subset_2_train", 
                                          "ARIMA"), 
                                color_types = c("ARIMA" = "grey", 
                                                "subset_2_gabriel" = "#00BFC4", 
                                                "subset_2_relative" = "#00B0F6", 
                                                "subset_2_soi" = "#9590FF", 
                                                "subset_2_delaunay" = "#E76BF3", 
                                                "subset_2_train" = "#FF62BC"))

m_2_knn_I <- plot_mase_II(mase_overview = mase_2_overview,
                          mase_name = "free", 
                          counties_subset = all_counties[1:9],
                          number_counties = 1, 
                          types = c("subset_2_dnn", 
                                    "subset_2_knn", 
                                    "subset_2_queen", 
                                    "subset_2_eco_hub", 
                                    "subset_2_complete", 
                                    "ARIMA"), 
                          color_types = c("ARIMA" = "grey", 
                                          "subset_2_knn" = "#F8766D", 
                                          "subset_2_dnn" = "#D89000", 
                                          "subset_2_complete" = "#A3A500", 
                                          "subset_2_queen" = "#39B600", 
                                          "subset_2_eco_hub" = "#00BF7D"))
m_2_knn_II <- plot_mase_II(mase_overview = mase_2_overview, 
                           mase_name = "free", 
                           counties_subset = all_counties[10:18],
                           number_counties = 2, 
                           types = c("subset_2_dnn", 
                                     "subset_2_knn", 
                                     "subset_2_queen", 
                                     "subset_2_eco_hub", 
                                     "subset_2_complete", 
                                     "ARIMA"), 
                           color_types = c("ARIMA" = "grey", 
                                           "subset_2_knn" = "#F8766D", 
                                           "subset_2_dnn" = "#D89000", 
                                           "subset_2_complete" = "#A3A500", 
                                           "subset_2_queen" = "#39B600", 
                                           "subset_2_eco_hub" = "#00BF7D"))
m_2_knn_III <- plot_mase_II(mase_overview = mase_2_overview,
                            mase_name = "free", 
                            counties_subset = all_counties[19:26],
                            number_counties = 3, 
                            types = c("subset_2_dnn", 
                                      "subset_2_knn", 
                                      "subset_2_queen", 
                                      "subset_2_eco_hub", 
                                      "subset_2_complete", 
                                      "ARIMA"), 
                            color_types = c("ARIMA" = "grey", 
                                            "subset_2_knn" = "#F8766D", 
                                            "subset_2_dnn" = "#D89000", 
                                            "subset_2_complete" = "#A3A500", 
                                            "subset_2_queen" = "#39B600", 
                                            "subset_2_eco_hub" = "#00BF7D"))


# Predicted vs.  fitted for unrestricted ----------------------------------
g_2_delaunay_I <- plot_predicted_vs_fitted_I(mase_overview = mase_2_overview, 
                                             mase_name = "free", 
                                             counties_subset = all_counties[1:9],
                                             number_counties = 1, 
                                             types = c("subset_2_delaunay", 
                                                       "subset_2_gabriel", 
                                                       "subset_2_relative", 
                                                       "subset_2_soi", 
                                                       "subset_2_train", 
                                                       "ARIMA"), 
                                             color_types = c("ARIMA" = "grey", 
                                                             "subset_2_gabriel" = "#00BFC4", 
                                                             "subset_2_relative" = "#00B0F6", 
                                                             "subset_2_soi" = "#9590FF", 
                                                             "subset_2_delaunay" = "#E76BF3", 
                                                             "subset_2_train" = "#FF62BC"))
g_2_delaunay_II <- plot_predicted_vs_fitted_I(mase_overview = mase_2_overview,
                                              mase_name = "free", 
                                              counties_subset = all_counties[10:18], 
                                              number_counties = 2, 
                                              types = c("subset_2_delaunay", 
                                                        "subset_2_gabriel", 
                                                        "subset_2_relative", 
                                                        "subset_2_soi", 
                                                        "subset_2_train", 
                                                        "ARIMA"), 
                                              color_types = c("ARIMA" = "grey", 
                                                              "subset_2_gabriel" = "#00BFC4", 
                                                              "subset_2_relative" = "#00B0F6", 
                                                              "subset_2_soi" = "#9590FF", 
                                                              "subset_2_delaunay" = "#E76BF3", 
                                                              "subset_2_train" = "#FF62BC"))
g_2_delaunay_III <- plot_predicted_vs_fitted_I(mase_overview = mase_2_overview, 
                                               mase_name = "free", 
                                               counties_subset = all_counties[19:26], 
                                               number_counties = 3, 
                                               types = c("subset_2_delaunay", 
                                                         "subset_2_gabriel", 
                                                         "subset_2_relative", 
                                                         "subset_2_soi", 
                                                         "subset_2_train", 
                                                         "ARIMA"), 
                                               color_types = c("ARIMA" = "grey", 
                                                               "subset_2_gabriel" = "#00BFC4", 
                                                               "subset_2_relative" = "#00B0F6", 
                                                               "subset_2_soi" = "#9590FF", 
                                                               "subset_2_delaunay" = "#E76BF3", 
                                                               "subset_2_train" = "#FF62BC"))

g_2_knn_I <- plot_predicted_vs_fitted_II(mase_overview = mase_2_overview,
                                         mase_name = "free",
                                         counties_subset = all_counties[1:9],
                                         number_counties = 1, 
                                         types = c("subset_2_dnn", 
                                                   "subset_2_knn", 
                                                   "subset_2_queen", 
                                                   "subset_2_eco_hub", 
                                                   "subset_2_complete", 
                                                   "ARIMA"), 
                                         color_types = c("ARIMA" = "grey", 
                                                         "subset_2_knn" = "#F8766D", 
                                                         "subset_2_dnn" = "#D89000", 
                                                         "subset_2_complete" = "#A3A500", 
                                                         "subset_2_queen" = "#39B600", 
                                                         "subset_2_eco_hub" = "#00BF7D"))
g_2_knn_II <- plot_predicted_vs_fitted_II(mase_overview = mase_2_overview,
                                          mase_name = "free",
                                          counties_subset = all_counties[10:18],
                                          number_counties = 2, 
                                          types = c("subset_2_dnn", 
                                                    "subset_2_knn", 
                                                    "subset_2_queen", 
                                                    "subset_2_eco_hub", 
                                                    "subset_2_complete", 
                                                    "ARIMA"), 
                                          color_types = c("ARIMA" = "grey", 
                                                          "subset_2_knn" = "#F8766D", 
                                                          "subset_2_dnn" = "#D89000", 
                                                          "subset_2_complete" = "#A3A500", 
                                                          "subset_2_queen" = "#39B600", 
                                                          "subset_2_eco_hub" = "#00BF7D"))
g_2_knn_III <- plot_predicted_vs_fitted_II(mase_overview = mase_2_overview,
                                           mase_name = "free",
                                           counties_subset = all_counties[19:26],
                                           number_counties = 3, 
                                           types = c("subset_2_dnn", 
                                                     "subset_2_knn", 
                                                     "subset_2_queen", 
                                                     "subset_2_eco_hub", 
                                                     "subset_2_complete", 
                                                     "ARIMA"), 
                                           color_types = c("ARIMA" = "grey", 
                                                           "subset_2_knn" = "#F8766D", 
                                                           "subset_2_dnn" = "#D89000", 
                                                           "subset_2_complete" = "#A3A500", 
                                                           "subset_2_queen" = "#39B600", 
                                                           "subset_2_eco_hub" = "#00BF7D"))

# Best GNAR for restrictions ----------------------------------------------
best_for_subset_all

best_model_restrictive <- fit_and_predict(alpha = 5, 
                                          beta = c(2, 1, 1, 1, 1), 
                                          net = covid_net_queen_gnar, 
                                          vts = datasets_list_coarse[[1]], 
                                          globalalpha = TRUE, 
                                          old = TRUE,
                                          forecast_window = 5, 
                                          return_model = TRUE)
# data set 2
best_model_free <- fit_and_predict(alpha = 5, 
                                   beta = c(1, 1, 1, 1, 0), 
                                   net = knn_21_gnar, 
                                   vts = datasets_list_coarse[[2]], 
                                   globalalpha = TRUE, 
                                   old = TRUE,
                                   forecast_window = 5, 
                                   return_model = TRUE)


mean_restrictive <- mase_1_queen %>% 
  dplyr::select(res, mase) %>% 
  colMeans()

mean_free <- mase_2_knn %>% 
  dplyr::select(res, mase) %>% 
  colMeans()


best_for_subset_all <- cbind(best_for_subset_all,
                             "AIC" = c(AIC(best_model_restrictive), 
                                       AIC(best_model_free)), 
                             rbind(mean_restrictive, 
                                   mean_free))

# for latex 
strCaption <- "Overview over the best performing model and network for 
restricted and unrestricted pandemic phases; average residual 
$\\Bar{\\varepsilon}$ and average (av.) MASE indicated for the predicted 5 
weeks at the end of the observed time period, 11.04.2021 - 09.05.2021 for the 
restricted dataset and 25.12.2022 - 22.01.2023 for the unrestricted dataset."
print(xtable(best_for_subset_all[, c(1, 4, 2, 3, 6, 7)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_pandemic_phases", 
             align = c("", "l", "|", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_for_subset_all[, c(1, 4, 2, 3, 6, 7)])),
                        command = c(paste("\\toprule \n",
                                          "Data subset & network & \\code{GNAR} model & 
                                          BIC  & $\\Bar{\\varepsilon}$ & av. MASE \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

paste0(best_for_subset_all$AIC %>% round(2), 
       collapse = ", ")



mean_arima_df <- data.frame(subset = c("Restricted", 
                              "Unrestricted"), 
                            net = rep("ARIMA", 2),
                            model = rep("", 2)) %>% 
  cbind(rbind(mase_arima_restrictive %>% 
               dplyr::select(-c(CountyName, true, type, time)) %>% 
               colMeans(), 
             mase_arima_free %>% 
               dplyr::select(-c(CountyName, true, type, time)) %>% 
               colMeans()), 
        mean_arima_BIC)

# for latex
print(xtable(mean_arima_df[, c(1, 2, 3, 7, 5, 6)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_pandemic_phases", 
             align = c("", "l", "|", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(mean_arima_df[, c(1, 2, 3, 7, 5, 6)])),
                        command = c(paste("\\toprule \n",
                                          "Data subset & network & \\code{GNAR} model & 
                                          BIC  & $\\Bar{\\varepsilon}$ & av. MASE \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

mean_arima_df$mean_AIC



# Residual analysis for pandemic phases -----------------------------------
residuals_restrictive <- check_and_plot_residuals_subset(model = best_model_restrictive, 
                                                         network_name = "restrictive_1_lag", 
                                                         alpha = 5, 
                                                         n_ahead = 5, 
                                                         data = datasets_list_coarse[[1]], 
                                                         counties = all_counties, 
                                                         dataset_name = "GNAR_pandemic_phases")

ks_1_significant <- mase_1_queen %>% 
  split(mase_1_queen$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value <= 0.025)
  }) %>% 
  list.cbind() %>% 
  table()

ks_1 <- mase_1_queen %>% 
  split(mase_1_queen$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value %>% round(digits = 3))
  }) %>% 
  list.rbind() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "CountyName")

residuals_free <- check_and_plot_residuals_subset(model = best_model_free, 
                                                  network_name = "free_1_lag", 
                                                  alpha = 5, 
                                                  n_ahead = 5, 
                                                  data = datasets_list_coarse[[2]], 
                                                  counties = all_counties, 
                                                  dataset_name = "GNAR_pandemic_phases")

ks_2_significant <- mase_2_knn %>% 
  split(mase_2_knn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value <= 0.025)
  }) %>% 
  list.cbind() %>% 
  table()

ks_2 <- mase_2_knn %>% 
  split(mase_2_knn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value %>% round(digits = 3))
  }) %>% 
  list.rbind() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "CountyName")

# MASE for each county ----------------------------------------------------

mean_restrictive_counties <- mase_1_queen %>% 
  dplyr::select(CountyName, res, mase) %>% 
  group_by(CountyName) %>% 
  summarise(res = mean(res), 
            mase = mean(mase))

overview_restrictive_counties <- left_join(mean_restrictive_counties, 
                                           ks_1, 
                                           by = "CountyName")

mean_free_counties <- mase_2_knn %>% 
  dplyr::select(CountyName, res, mase) %>% 
  group_by(CountyName) %>% 
  summarise(res = mean(res), 
            mase = mean(mase))

overview_free_counties <- left_join(mean_free_counties, 
                                           ks_2, 
                                           by = "CountyName")


overview_counties <- left_join(overview_restrictive_counties, overview_free_counties, 
                               by = "CountyName")

# for latex 
strCaption <- "Average residual $\\Bar{\\varepsilon}$ and average (av.) MASE value 
as well as Kolmogorov-Smirnov p-value (p) for each county for restricted and unrestricted 
pandemic phase"
print(xtable(overview_counties,
             digits=2,
             caption=strCaption,
             label="tab:mase_res_p_counties", 
             align = c("", "l", "|", "r", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(overview_counties)),
                        command = c(paste("\\toprule \n",
                                          "County & $\\Bar{\\varepsilon}$ & 
                                          av. MASE & p & 
                                          $\\Bar{\\varepsilon}$ & av. MASE & 
                                          p \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Change in coefficients --------------------------------------------------
best_for_subset_all

parameter_development_phases(net_list = list(covid_net_queen_gnar, 
                                             knn_21_gnar), 
                             alpha_vector = c(5, 5),
                             beta_list = list(c(2, 1, 1, 1, 1), 
                                              c(1, 1, 1, 1, 0)),
                             county_index = county_index_knn,
                             globalalpha = TRUE,
                             old = TRUE)

param_knn <- parameter_development_phases(data_list = datasets_list_coarse , 
                                          net_list = list(knn_21_gnar, 
                                                          knn_21_gnar), 
                                          alpha = c(5, 5), 
                                          beta = list(c(1, 1, 1, 1, 0),
                                                      c(1, 1, 1, 1, 0)), 
                                          globalalpha = TRUE, 
                                          old = TRUE, 
                                          name = "phases_same_model")
