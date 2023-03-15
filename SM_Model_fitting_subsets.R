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

# Load data ---------------------------------------------------------------
load(file = "Data/RObjects/data_subsets_regulations.RData")

# Best performing model for every subset ----------------------------------
# compute GNAR models for each data subsets and select the best performing one 
# based on the BIC 

# Queen 
best_for_subset_queen <- fit_and_predict_for_restrictions(net = covid_net_queen_gnar)
best_for_subset_queen$network <- "Queen"

# Economic hub
best_for_subset_eco_hub <- fit_and_predict_for_restrictions(net = covid_net_eco_hubs_gnar, 
                                                            numeric_vertices = TRUE, 
                                                            county_index = county_index_eco_hubs, 
                                                            upper_limit = 4)
best_for_subset_eco_hub$network <- "Eco. hub"

# Railway-based 
best_for_subset_train <- fit_and_predict_for_restrictions(net = covid_net_train_gnar, 
                                                          numeric_vertices = TRUE, 
                                                          county_index = county_index_train)
best_for_subset_train$network <- "Train"

# Delaunay triangulation 
best_for_subset_delaunay <- fit_and_predict_for_restrictions(net = covid_net_delaunay_gnar)
best_for_subset_delaunay$network <- "Delaunay"

# Gabriel 
best_for_subset_gabriel <- fit_and_predict_for_restrictions(net = covid_net_gabriel_gnar)
best_for_subset_gabriel$network <- "Gabriel"

# Relative neighbourhood
best_for_subset_relative <- fit_and_predict_for_restrictions(net = covid_net_relative_gnar)
best_for_subset_relative$network <- "Relative"

# SOI 
best_for_subset_soi <- fit_and_predict_for_restrictions(net = covid_net_soi_gnar)
best_for_subset_soi$network <- "SOI"

# Complete
best_for_subset_complete <- fit_and_predict_for_restrictions(net = complete_net_gnar, 
                                                             upper_limit = 1)
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
                                          upper_limit = upper_limit_knn)
  
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
                                          upper_limit = upper_limit_dnn)
  
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



# Construct networks ------------------------------------------------------
# construct the networks for the best performing GNAR models 

# DNN d = 300
dnn_300 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 300,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_300_igraph<- neighborsDataFrame(nb = dnn_300) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_300_gnar <- dnn_300_igraph %>% 
  igraphtoGNAR()

# DNN d = 175
dnn_175 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 175,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_175_igraph<- neighborsDataFrame(nb = dnn_175) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_175_gnar <- dnn_175_igraph %>% 
  igraphtoGNAR()

# DNN d = 225 
dnn_225 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 225,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_225_igraph<- neighborsDataFrame(nb = dnn_225) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_225_gnar <- dnn_225_igraph %>% 
  igraphtoGNAR()

# DNN d = 100 
dnn_100 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 100,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_100_igraph<- neighborsDataFrame(nb = dnn_100) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_100_gnar <- dnn_100_igraph %>% 
  igraphtoGNAR()

# DNN d = 275
dnn_275 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 275,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_275_igraph<- neighborsDataFrame(nb = dnn_275) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_275_gnar <- dnn_275_igraph %>% 
  igraphtoGNAR()

# DNN d = 150 
dnn_150 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 150,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_150_igraph<- neighborsDataFrame(nb = dnn_150) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_150_gnar <- dnn_150_igraph %>% 
  igraphtoGNAR()

# DNN d = 250
dnn_250 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 250,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_250_igraph<- neighborsDataFrame(nb = dnn_250) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_250_gnar <- dnn_250_igraph %>% 
  igraphtoGNAR()


# KNN k = 1
knn_1 <- knearneigh(x = coord_urbanisation, 
                    k = 1, 
                    longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_1_igraph<- neighborsDataFrame(nb = knn_1) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_1_gnar <- knn_1_igraph %>% 
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

# KNN k = 5
knn_17 <- knearneigh(x = coord_urbanisation, 
                     k = 17, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_17_igraph<- neighborsDataFrame(nb = knn_17) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_17_gnar <- knn_17_igraph %>% 
  igraphtoGNAR()

# KNN k = 9
knn_9 <- knearneigh(x = coord_urbanisation, 
                    k = 9, 
                    longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_9_igraph<- neighborsDataFrame(nb = knn_9) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_9_gnar <- knn_9_igraph %>% 
  igraphtoGNAR()

# KNN k = 11
knn_11 <- knearneigh(x = coord_urbanisation, 
                     k = 11, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_11_igraph<- neighborsDataFrame(nb = knn_11) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_11_gnar <- knn_11_igraph %>% 
  igraphtoGNAR()


# KNN k = 15
knn_15 <- knearneigh(x = coord_urbanisation, 
                     k = 15, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_15_igraph<- neighborsDataFrame(nb = knn_15) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_15_gnar <- knn_15_igraph %>% 
  igraphtoGNAR()

# KNN k = 25
knn_25 <- knearneigh(x = coord_urbanisation, 
                     k = 25, 
                     longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_25_igraph<- neighborsDataFrame(nb = knn_25) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_25_gnar <- knn_25_igraph %>% 
  igraphtoGNAR()


# Compare networks --------------------------------------------------------
# compute network characteristics

# DNN 175
char_dnn_175 <- network_characteristics(dnn_175_igraph, 
                                        network_name = "dnn_175")

# DNN 225
char_dnn_225 <- network_characteristics(dnn_225_igraph, 
                                        network_name = "dnn_225")

# DNN 300
char_dnn_300 <- network_characteristics(dnn_300_igraph, 
                                        network_name = "dnn_300")

# KNN 9
char_knn_9 <- network_characteristics(knn_9_igraph, 
                                      network_name = "knn_9")

# Economic hub 
char_eco_hubs <- network_characteristics(covid_net_eco_hubs_igraph, 
                                         network_name = "eco_hub")

# compare network characteristics 
cbind(char_dnn_300[, 2], 
      char_dnn_175, 
      char_knn_9[, 2], 
      char_dnn_225[, 2], 
      char_eco_hubs[, 2])



# Best GNAR model for each data subset ------------------------------------
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

# fit models 
# data set 1
best_model_1 <- fit_and_predict(alpha = 5, 
                                beta = c(1, 0, 0, 0, 0), 
                                net = dnn_300_gnar, 
                                vts = datasets_list[[1]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 2
best_model_2 <- fit_and_predict(alpha = 5, 
                                beta = c(1, 1, 1, 1, 1), 
                                net = dnn_175_gnar, 
                                vts = datasets_list[[2]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 3
best_model_3 <- fit_and_predict(alpha = 5, 
                                beta = c(4, 0, 0, 0, 0), 
                                net = covid_net_queen_gnar, 
                                vts = datasets_list[[3]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 4
best_model_4 <- fit_and_predict(alpha = 3, 
                                beta = c(1, 1, 1), 
                                net = dnn_250_gnar, 
                                vts = datasets_list[[4]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 5
best_model_5 <- fit_and_predict(alpha = 5, 
                                beta = c(4, 2, 2, 2, 0), 
                                net = covid_net_gabriel_gnar, 
                                vts = datasets_list[[5]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)

best_for_subset_all$AIC <- c(AIC(best_model_1), 
                             AIC(best_model_2), 
                             AIC(best_model_3), 
                             AIC(best_model_4), 
                             AIC(best_model_5))

# Data formatting ----------------------------------------------------------
# format data subsets as data frame with time column 
data_1 <- datasets_list[[1]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[1]]) %>% as.Date())

data_2 <- datasets_list[[2]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[2]])  %>% as.Date())

data_3 <- datasets_list[[3]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[3]])  %>% as.Date())

data_4 <- datasets_list[[4]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[4]])  %>% as.Date())

data_5 <- datasets_list[[5]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[5]])  %>% as.Date())



# MASE for data subset 1 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 1
best_for_subset %>% filter(data_subset == 1)
best_subset_knn_dnn %>% filter(data_subset == 1)

# Queen
mod_1_queen <- fit_and_predict(alpha = 5, 
                               beta = c(3, 2, 2, 1, 1), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_1_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 1, 1, 0, 0), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_1_train <- fit_and_predict(alpha = 5, 
                               beta = c(1, 0, 0, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_1_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(2, 2, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_1_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(3, 1, 1, 0, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_1_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(4, 0, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_1_soi <- fit_and_predict(alpha = 5, 
                             beta = c(4, 1, 1, 1, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_1_knn <- fit_and_predict(alpha = 5, 
                             beta = c(4, 1, 1, 0, 0), 
                             net = knn_7_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_1_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 0, 0, 0, 0), 
                             net = dnn_300_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_1_complete <- fit_and_predict(alpha = 4, 
                                  beta = c(1, 1, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_1_queen <- compute_MASE(model = mod_1_queen, 
                             network_name = "subset_1_queen", 
                             n_ahead = 5, 
                             data_df = data_1)

mase_1_eco_hub <- compute_MASE(model = mod_1_eco_hub, 
                               network_name = "subset_1_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_1)

mase_1_train <- compute_MASE(model = mod_1_train, 
                             network_name = "subset_1_train", 
                             n_ahead = 5, 
                             data_df = data_1)

mase_1_delaunay <- compute_MASE(model = mod_1_delaunay, 
                                network_name = "subset_1_delaunay", 
                                n_ahead = 5, 
                                data_df = data_1)

mase_1_gabriel <- compute_MASE(model = mod_1_gabriel, 
                               network_name = "subset_1_gabriel", 
                               n_ahead = 5, 
                               data_df = data_1)

mase_1_relative <- compute_MASE(model = mod_1_relative, 
                                network_name = "subset_1_relative", 
                                n_ahead = 5, 
                                data_df = data_1)

mase_1_soi <- compute_MASE(model = mod_1_soi, 
                           network_name = "subset_1_soi", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_knn <- compute_MASE(model = mod_1_knn, 
                           network_name = "subset_1_knn", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_dnn <- compute_MASE(model = mod_1_dnn, 
                           network_name = "subset_1_dnn", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_complete <- compute_MASE(model = mod_1_complete, 
                                network_name = "subset_1_complete", 
                                n_ahead = 5, 
                                data_df = data_1)


mase_1_overview <- rbind(mase_1_queen, 
                         mase_1_eco_hub, 
                         mase_1_train, 
                         mase_1_delaunay, 
                         mase_1_gabriel, 
                         mase_1_relative, 
                         mase_1_soi, 
                         mase_1_knn, 
                         mase_1_dnn, 
                         mase_1_complete)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network
ggplot(mase_1_overview %>% filter(type %in% c("subset_1_delaunay", 
                                              "subset_1_gabriel", 
                                              "subset_1_relative", 
                                              "subset_1_soi", 
                                              "subset_1_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_1_gabriel" = "#00BFC4", 
                                "subset_1_relative" = "#00B0F6", 
                                "subset_1_soi" = "#9590FF", 
                                "subset_1_delaunay" = "#E76BF3", 
                                "subset_1_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_delaunay_etc_subset_1.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_1_overview %>% filter(type %in% c("subset_1_dnn", 
                                              "subset_1_knn", 
                                              "subset_1_queen", 
                                              "subset_1_eco_hub", 
                                              "subset_1_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_1_knn" = "#F8766D", 
                                "subset_1_dnn" = "#D89000", 
                                "subset_1_complete" = "#A3A500", 
                                "subset_1_queen" = "#39B600", 
                                "subset_1_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_knn_etc_subset_1.pdf", 
       width = 26, height = 13, units = "cm")





# MASE for data subset 2 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 2
best_for_subset %>% filter(data_subset == 2)
best_subset_knn_dnn %>% filter(data_subset == 2)

# Queen
mod_2_queen <- fit_and_predict(alpha = 5, 
                               beta = c(2, 2, 1, 1, 0), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[2]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_2_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 2, 2, 1, 1), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[2]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_2_train <- fit_and_predict(alpha = 5, 
                               beta = c(5, 1, 1, 1, 1), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[2]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_2_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(3, 1, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_2_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(2, 2, 2, 2, 2), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[2]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_2_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(3, 0, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_2_soi <- fit_and_predict(alpha = 5, 
                             beta = c(2, 1, 1, 1, 1), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_2_knn <- fit_and_predict(alpha = 5, 
                             beta = c(2, 2, 2, 2, 1), 
                             net = knn_7_gnar, 
                             vts = datasets_list[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_2_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 1), 
                             net = dnn_175_gnar, 
                             vts = datasets_list[[2]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_2_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 0, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[2]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_2_queen <- compute_MASE(model = mod_2_queen, 
                             network_name = "subset_2_queen", 
                             n_ahead = 5, 
                             data_df = data_2)

mase_2_eco_hub <- compute_MASE(model = mod_2_eco_hub, 
                               network_name = "subset_2_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_2)

mase_2_train <- compute_MASE(model = mod_2_train, 
                             network_name = "subset_2_train", 
                             n_ahead = 5, 
                             data_df = data_2)

mase_2_delaunay <- compute_MASE(model = mod_2_delaunay, 
                                network_name = "subset_2_delaunay", 
                                n_ahead = 5, 
                                data_df = data_2)

mase_2_gabriel <- compute_MASE(model = mod_2_gabriel, 
                               network_name = "subset_2_gabriel", 
                               n_ahead = 5, 
                               data_df = data_2)

mase_2_relative <- compute_MASE(model = mod_2_relative, 
                                network_name = "subset_2_relative", 
                                n_ahead = 5, 
                                data_df = data_2)

mase_2_soi <- compute_MASE(model = mod_2_soi, 
                           network_name = "subset_2_soi", 
                           n_ahead = 5, 
                           data_df = data_2)

mase_2_knn <- compute_MASE(model = mod_2_knn, 
                           network_name = "subset_2_knn", 
                           n_ahead = 5, 
                           data_df = data_2)

mase_2_dnn <- compute_MASE(model = mod_2_dnn, 
                           network_name = "subset_2_dnn", 
                           n_ahead = 5, 
                           data_df = data_2)

mase_2_complete <- compute_MASE(model = mod_2_complete, 
                                network_name = "subset_2_complete", 
                                n_ahead = 5, 
                                data_df = data_2)


mase_2_overview <- rbind(mase_2_queen, 
                         mase_2_eco_hub, 
                         mase_2_train, 
                         mase_2_delaunay, 
                         mase_2_gabriel, 
                         mase_2_relative, 
                         mase_2_soi, 
                         mase_2_knn, 
                         mase_2_dnn, 
                         mase_2_complete)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network
ggplot(mase_2_overview %>% filter(type %in% c("subset_2_delaunay", 
                                              "subset_2_gabriel", 
                                              "subset_2_relative", 
                                              "subset_2_soi", 
                                              "subset_2_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_2_gabriel" = "#00BFC4", 
                                "subset_2_relative" = "#00B0F6", 
                                "subset_2_soi" = "#9590FF", 
                                "subset_2_delaunay" = "#E76BF3", 
                                "subset_2_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_delaunay_etc_subset_2.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_2_overview %>% filter(type %in% c("subset_2_dnn", 
                                              "subset_2_knn", 
                                              "subset_2_queen", 
                                              "subset_2_eco_hub", 
                                              "subset_2_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_2_knn" = "#F8766D", 
                                "subset_2_dnn" = "#D89000", 
                                "subset_2_complete" = "#A3A500", 
                                "subset_2_queen" = "#39B600", 
                                "subset_2_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_knn_etc_subset_2.pdf", 
       width = 26, height = 13, units = "cm")




# MASE for data subset 3 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 3
best_for_subset %>% filter(data_subset == 3)
best_subset_knn_dnn %>% filter(data_subset == 3)

# Queen
mod_3_queen <- fit_and_predict(alpha = 5, 
                               beta = c(4, 0, 0, 0, 0), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[3]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_3_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(3, 1, 0, 0, 0), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[3]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_3_train <- fit_and_predict(alpha = 5, 
                               beta = c(2, 2, 2, 2, 1), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[3]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_3_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(5, 1, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[3]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_3_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(2, 2, 2, 2, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[3]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_3_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(4, 1, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[3]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_3_soi <- fit_and_predict(alpha = 5, 
                             beta = c(2, 1, 0, 0, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[3]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_3_knn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 0, 0, 0, 0), 
                             net = knn_15_gnar, 
                             vts = datasets_list[[3]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_3_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(2, 2, 2, 2, 1), 
                             net = dnn_100_gnar, 
                             vts = datasets_list[[3]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_3_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 1, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[3]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_3_queen <- compute_MASE(model = mod_3_queen, 
                             network_name = "subset_3_queen", 
                             n_ahead = 5, 
                             data_df = data_3)

mase_3_eco_hub <- compute_MASE(model = mod_3_eco_hub, 
                               network_name = "subset_3_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_3)

mase_3_train <- compute_MASE(model = mod_3_train, 
                             network_name = "subset_3_train", 
                             n_ahead = 5, 
                             data_df = data_3)

mase_3_delaunay <- compute_MASE(model = mod_3_delaunay, 
                                network_name = "subset_3_delaunay", 
                                n_ahead = 5, 
                                data_df = data_3)

mase_3_gabriel <- compute_MASE(model = mod_3_gabriel, 
                               network_name = "subset_3_gabriel", 
                               n_ahead = 5, 
                               data_df = data_3)

mase_3_relative <- compute_MASE(model = mod_3_relative, 
                                network_name = "subset_3_relative", 
                                n_ahead = 5, 
                                data_df = data_3)

mase_3_soi <- compute_MASE(model = mod_3_soi, 
                           network_name = "subset_3_soi", 
                           n_ahead = 5, 
                           data_df = data_3)

mase_3_knn <- compute_MASE(model = mod_3_knn, 
                           network_name = "subset_3_knn", 
                           n_ahead = 5, 
                           data_df = data_3)

mase_3_dnn <- compute_MASE(model = mod_3_dnn, 
                           network_name = "subset_3_dnn", 
                           n_ahead = 5, 
                           data_df = data_3)

mase_3_complete <- compute_MASE(model = mod_3_complete, 
                                network_name = "subset_3_complete", 
                                n_ahead = 5, 
                                data_df = data_3)


mase_3_overview <- rbind(mase_3_queen, 
                         mase_3_eco_hub, 
                         mase_3_train, 
                         mase_3_delaunay, 
                         mase_3_gabriel, 
                         mase_3_relative, 
                         mase_3_soi, 
                         mase_3_knn, 
                         mase_3_dnn, 
                         mase_3_complete)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network
ggplot(mase_3_overview %>% filter(type %in% c("subset_3_delaunay", 
                                              "subset_3_gabriel", 
                                              "subset_3_relative", 
                                              "subset_3_soi", 
                                              "subset_3_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_3_gabriel" = "#00BFC4", 
                                "subset_3_relative" = "#00B0F6", 
                                "subset_3_soi" = "#9590FF", 
                                "subset_3_delaunay" = "#E76BF3", 
                                "subset_3_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_delaunay_etc_subset_3.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_3_overview %>% filter(type %in% c("subset_3_dnn", 
                                              "subset_3_knn", 
                                              "subset_3_queen", 
                                              "subset_3_eco_hub", 
                                              "subset_3_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_3_knn" = "#F8766D", 
                                "subset_3_dnn" = "#D89000", 
                                "subset_3_complete" = "#A3A500", 
                                "subset_3_queen" = "#39B600", 
                                "subset_3_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_knn_etc_subset_3.pdf", 
       width = 26, height = 13, units = "cm")




# MASE for data subset 4 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 4
best_for_subset %>% filter(data_subset == 4)
best_subset_knn_dnn %>% filter(data_subset == 4)

# Queen
mod_4_queen <- fit_and_predict(alpha = 2, 
                               beta = c(2, 2), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[4]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)


# Eco hub 
mod_4_eco_hub <- fit_and_predict(alpha = 2, 
                                 beta = c(1, 1), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[4]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_4_train <- fit_and_predict(alpha = 5, 
                               beta = c(5, 0, 0, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[4]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_4_delaunay <- fit_and_predict(alpha = 3, 
                                  beta = c(5, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_4_gabriel <- fit_and_predict(alpha = 2, 
                                 beta = c(2, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[4]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_4_relative <- fit_and_predict(alpha = 3, 
                                  beta = c(5, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_4_soi <- fit_and_predict(alpha = 3, 
                             beta = c(4, 0, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_4_knn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 1), 
                             net = knn_25_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_4_dnn <- fit_and_predict(alpha = 3, 
                             beta = c(1, 1, 1), 
                             net = dnn_250_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_4_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# compute MASE for best performing GNAR model for each network 
mase_4_queen <- compute_MASE(model = mod_4_queen, 
                             network_name = "subset_4_queen", 
                             n_ahead = 5, 
                             data_df = data_4)

mase_4_eco_hub <- compute_MASE(model = mod_4_eco_hub, 
                               network_name = "subset_4_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_4)

mase_4_train <- compute_MASE(model = mod_4_train, 
                             network_name = "subset_4_train", 
                             n_ahead = 5, 
                             data_df = data_4)

mase_4_delaunay <- compute_MASE(model = mod_4_delaunay, 
                                network_name = "subset_4_delaunay", 
                                n_ahead = 5, 
                                data_df = data_4)

mase_4_gabriel <- compute_MASE(model = mod_4_gabriel, 
                               network_name = "subset_4_gabriel", 
                               n_ahead = 5, 
                               data_df = data_4)

mase_4_relative <- compute_MASE(model = mod_4_relative, 
                                network_name = "subset_4_relative", 
                                n_ahead = 5, 
                                data_df = data_4)

mase_4_soi <- compute_MASE(model = mod_4_soi, 
                           network_name = "subset_4_soi", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_knn <- compute_MASE(model = mod_4_knn, 
                           network_name = "subset_4_knn", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_dnn <- compute_MASE(model = mod_4_dnn, 
                           network_name = "subset_4_dnn", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_complete <- compute_MASE(model = mod_4_complete, 
                                network_name = "subset_4_complete", 
                                n_ahead = 5, 
                                data_df = data_4)


mase_4_overview <- rbind(mase_4_queen, 
                         mase_4_eco_hub, 
                         mase_4_train, 
                         mase_4_delaunay, 
                         mase_4_gabriel, 
                         mase_4_relative, 
                         mase_4_soi, 
                         mase_4_knn, 
                         mase_4_dnn, 
                         mase_4_complete)

# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network 
ggplot(mase_4_overview %>% filter(type %in% c("subset_4_delaunay", 
                                              "subset_4_gabriel", 
                                              "subset_4_relative", 
                                              "subset_4_soi", 
                                              "subset_4_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_4_gabriel" = "#00BFC4", 
                                "subset_4_relative" = "#00B0F6", 
                                "subset_4_soi" = "#9590FF", 
                                "subset_4_delaunay" = "#E76BF3", 
                                "subset_4_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_delaunay_etc_subset_4.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_4_overview %>% filter(type %in% c("subset_4_dnn", 
                                              "subset_4_knn", 
                                              "subset_4_queen", 
                                              "subset_4_eco_hub", 
                                              "subset_4_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_4_knn" = "#F8766D", 
                                "subset_4_dnn" = "#D89000", 
                                "subset_4_complete" = "#A3A500", 
                                "subset_4_queen" = "#39B600", 
                                "subset_4_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_knn_etc_subset_4.pdf", 
       width = 26, height = 13, units = "cm")



# MASE for data subset 5 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 5
best_for_subset %>% filter(data_subset == 5)
best_subset_knn_dnn %>% filter(data_subset == 5)

# Queen
mod_5_queen <- fit_and_predict(alpha = 5, 
                               beta = c(4, 0, 0, 0, 0), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[5]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_5_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 0, 0, 0, 0),
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[5]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_5_train <- fit_and_predict(alpha = 5, 
                               beta = c(3, 1, 1, 1, 0),
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[5]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_5_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(4, 2, 2, 2, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[5]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_5_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 2, 2, 2, 0),
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[5]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_5_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(5, 0, 0, 0, 0),
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[5]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_5_soi <- fit_and_predict(alpha = 5, 
                             beta =  c(5, 1, 1, 1, 0),
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[5]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_5_knn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 0, 0),
                             net = knn_17_gnar, 
                             vts = datasets_list[[5]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_5_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 1, 0),
                             net = dnn_300_gnar, 
                             vts = datasets_list[[5]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_5_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 0, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[5]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_5_queen <- compute_MASE(model = mod_5_queen, 
                             network_name = "subset_5_queen", 
                             n_ahead = 5, 
                             data_df = data_5)

mase_5_eco_hub <- compute_MASE(model = mod_5_eco_hub, 
                               network_name = "subset_5_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_5)

mase_5_train <- compute_MASE(model = mod_5_train, 
                             network_name = "subset_5_train", 
                             n_ahead = 5, 
                             data_df = data_5)

mase_5_delaunay <- compute_MASE(model = mod_5_delaunay, 
                                network_name = "subset_5_delaunay", 
                                n_ahead = 5, 
                                data_df = data_5)

mase_5_gabriel <- compute_MASE(model = mod_5_gabriel, 
                               network_name = "subset_5_gabriel", 
                               n_ahead = 5, 
                               data_df = data_5)

mase_5_relative <- compute_MASE(model = mod_5_relative, 
                                network_name = "subset_5_relative", 
                                n_ahead = 5, 
                                data_df = data_5)

mase_5_soi <- compute_MASE(model = mod_5_soi, 
                           network_name = "subset_5_soi", 
                           n_ahead = 5, 
                           data_df = data_5)

mase_5_knn <- compute_MASE(model = mod_5_knn, 
                           network_name = "subset_5_knn", 
                           n_ahead = 5, 
                           data_df = data_5)

mase_5_dnn <- compute_MASE(model = mod_5_dnn, 
                           network_name = "subset_5_dnn", 
                           n_ahead = 5, 
                           data_df = data_5)

mase_5_complete <- compute_MASE(model = mod_5_complete, 
                                network_name = "subset_5_complete", 
                                n_ahead = 5, 
                                data_df = data_5)


mase_5_overview <- rbind(mase_5_queen, 
                         mase_5_eco_hub, 
                         mase_5_train, 
                         mase_5_delaunay, 
                         mase_5_gabriel, 
                         mase_5_relative, 
                         mase_5_soi, 
                         mase_5_knn, 
                         mase_5_dnn, 
                         mase_5_complete)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network
ggplot(mase_5_overview %>% filter(type %in% c("subset_5_delaunay", 
                                              "subset_5_gabriel", 
                                              "subset_5_relative", 
                                              "subset_5_soi", 
                                              "subset_5_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_5_gabriel" = "#00BFC4", 
                                "subset_5_relative" = "#00B0F6", 
                                "subset_5_soi" = "#9590FF", 
                                "subset_5_delaunay" = "#E76BF3", 
                                "subset_5_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_delaunay_etc_subset_5.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_5_overview %>% filter(type %in% c("subset_5_dnn", 
                                              "subset_5_knn", 
                                              "subset_5_queen", 
                                              "subset_5_eco_hub", 
                                              "subset_5_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_5_knn" = "#F8766D", 
                                "subset_5_dnn" = "#D89000", 
                                "subset_5_complete" = "#A3A500", 
                                "subset_5_queen" = "#39B600", 
                                "subset_5_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("Figures/GNAR_pandemic_phases/mase_knn_etc_subset_5.pdf", 
       width = 26, height = 13, units = "cm")



# MASE - Residual analysis ------------------------------------------------
best_for_subset_all
mean_1 <- mase_1_dnn %>% 
  dplyr::select(res, mase) %>% 
  tail(5) %>% 
  colMeans()

mean_2 <- mase_2_dnn %>% 
  dplyr::select(res, mase) %>% 
  tail(5) %>% 
  colMeans()

mean_3 <- mase_3_queen %>% 
  dplyr::select(res, mase) %>% 
  tail(5) %>% 
  colMeans()

mean_4 <- mase_4_dnn %>% 
  dplyr::select(res, mase) %>% 
  tail(5) %>% 
  colMeans()

mean_5 <- mase_5_gabriel %>% 
  dplyr::select(res, mase) %>% 
  tail(5) %>% 
  colMeans()

best_for_subset_all <- cbind(best_for_subset_all, 
                             as.data.frame(rbind(mean_1, 
                                                 mean_2, 
                                                 mean_3, 
                                                 mean_4, 
                                                 mean_5)))

# for latex 
strCaption <- "Overview over the best performing model and network for every 
COVID-19 data subset"
print(xtable(best_for_subset_all[, c(1, 4, 2, 3, 6, 7)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_datasets", 
             align = c("", "l", "|", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_for_subset_all[, c(1, 4, 2, 3, 6, 7)])),
                        command = c(paste("\\toprule \n",
                                          "Data subset & network & \\code{GNAR} model & 
                                          BIC  & mean residual & mean MASE \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

paste0(best_for_subset_all$AIC %>% round(2), 
       collapse = ", ")


# Residuals analysis ------------------------------------------------------
best_for_subset_all

# dataset 1
residuals_subset_1 <- check_and_plot_residuals_subset(model = mod_1_dnn, 
                                                      network_name = "subset_1_dnn", 
                                                      alpha = 5, 
                                                      n_ahead = 5, 
                                                      data = datasets_list[[1]])

# dataset 2
residuals_subset_2 <- check_and_plot_residuals_subset(model = mod_2_dnn, 
                                                      network_name = "subset_2_dnn", 
                                                      alpha = 5, 
                                                      n_ahead = 5, 
                                                      data = datasets_list[[2]])

# dataset 3
residuals_subset_3 <- check_and_plot_residuals_subset(model = mod_3_gabriel, 
                                                      network_name = "subset_3_knn", 
                                                      alpha = 5, 
                                                      n_ahead = 5, 
                                                      data = datasets_list[[3]])

# dataset 4
residuals_subset_4 <- check_and_plot_residuals_subset(model = mod_4_dnn, 
                                                      network_name = "subset_4_dnn", 
                                                      alpha = 5, 
                                                      n_ahead = 5, 
                                                      data = datasets_list[[4]])

# dataset 5
residuals_subset_5 <- check_and_plot_residuals_subset(model = mod_5_gabriel, 
                                                      network_name = "subset_5_eco_hub", 
                                                      alpha = 5, 
                                                      n_ahead = 5, 
                                                      data = datasets_list[[5]])


ks_1 <- mase_1_dnn %>% 
  split(mase_1_dnn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value)
  }) %>% 
  list.cbind()

ks_2 <- mase_2_dnn %>% 
  split(mase_2_dnn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value)
  }) %>% 
  list.cbind()

ks_3 <- mase_3_knn %>% 
  split(mase_3_knn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value)
  }) %>% 
  list.cbind()

ks_4 <- mase_4_dnn %>% 
  split(mase_4_dnn$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value)
  }) %>% 
  list.cbind()

ks_5 <- mase_5_gabriel %>% 
  split(mase_5_gabriel$CountyName) %>% 
  lapply(FUN = function(i) {
    return(ks.test(i$res, "pnorm")$p.value)
  }) %>% 
  list.cbind()



a  <- rbind(ks_1, 
            ks_2, 
            ks_3, 
            ks_4, 
            ks_5) %>% 
  round(4)



# Change in coefficients --------------------------------------------------
best_for_subset_all

parameter_development_subsets(net_list = list(dnn_300_gnar, 
                                              dnn_175_gnar, 
                                              covid_net_queen_gnar, 
                                              dnn_250_gnar, 
                                              covid_net_gabriel_gnar), 
                              alpha_vector = c(5, 5, 5, 3, 5),
                              beta_list = list(c(1, 0, 0, 0, 0), 
                                               c(1, 1, 1, 1, 1),
                                               c(4, 0, 0, 0, 0), 
                                               c(1, 1, 1), 
                                               c(4, 2, 2, 2, 0)),
                              county_index = county_index_dnn,
                              globalalpha = TRUE,
                              old = TRUE)


