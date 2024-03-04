### MODELLING COVID-19 INCIDENCE IN IRISH COUNTIES IN GNAR MODELS WITH 
### DIFFERENTLY CONSTRUCTED NETWORKS 

rm(list=ls())

# set seed to guarantee reproducibility 
set.seed(1234)

# Load libraries ----------------------------------------------------------
library(readr)
library(igraph)
library(GNAR)
library(MASS) # for box cox 
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 
library(geosphere) # for long / lat distance calculation
library(ade4) # igraph to neighbourhood list object
library(Hmisc) # for weighted variance 
library(Metrics) # for MASE computation 
library(rlist) # for easy concatenation of lists
library(ape)

# for MAPS: 
library(sp)
library(raster) 
library(leaflet)
library(mapview) # for saving

# for neighbourhood construction 
library(spdep) 

# for sphere of influence network 
# library(rgeos)
library(dbscan)

# for visualising the network
library(sf)
library(lwgeom)

# to transform to shapefile 
library(raster)

# for pairwise distances
library(flexclust)

# for latex expressions in graph titles 
library(latex2exp)

# for non-overlapping labels in base plot 
library(basicPlotteR)

# for arima models
library(forecast)


# load vectors and GNAR objects 
load(file = "Data/RObjects/GNAR.RData")
load(file = "Data/RObjects/igraph.RData")
load(file = "Data/RObjects/population_weight.RData")
load(file = "Data/RObjects/distance_urbanisation.RData")
load(file = "Data/RObjects/county_index.RData")
load(file = "Data/RObjects/coord_urbanisation.RData")
load(file = "Data/RObjects/urbanisation_factor.RData")

source("functions_paper.R")

# Turn off warnings
options(warn = -1)

# Set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load pre-processed data -------------------------------------------------
COVID_weekly_data <- read_csv(file = "Data/ireland_covid_weekly.csv", 
                       show_col_types = FALSE)

# correct format for GNAR models
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# transform to data frame with time column for plotting 
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 

# Weighting ---------------------------------------------------------------
# initialise data frame to store constructed weights and differences in 
# 1-lag COVID-19 ID
diff_all_counties <- data.frame("time" = NA, 
                                "CountyName" = NA, 
                                "cases" = NA, 
                                "dist" = NA,  
                                "pop_product" = NA, 
                                "pop_diff" = NA)
counties <- COVID_weekly_data$CountyName %>% unique()

for (county in counties) {
  # construct data frame for current county 
  cases_county <- covid_cases[, county]
  
  # compute difference in 1-lag COVID-19 ID to all the other counties 
  diff_cases <- covid_cases[, counties != county] - cases_county
  
  # format as a data frame with time column 
  diff_cases <- diff_cases %>% as.data.frame()
  diff_cases$time <- diff_cases %>% rownames()
  
  diff_cases_df <- diff_cases %>% gather("CountyName", "cases", -time)
  
  diff_county_names <- diff_cases %>% colnames()
  
  # filter row corresponding to the county in question from the
  # Great Circle distance data frame 
  dist_counties <- dist_urbanisation[county, 
                                     diff_county_names[diff_county_names != "time"]]
  
  # extract the population size for the county in question 
  county_pop <- population_weight[population_weight$CountyName == county, ]$weight
  
  
  
  dist_df <- data.frame("dist" = as.numeric(dist_counties[1, ]) / 1000, 
                        "CountyName" = dist_counties %>% colnames(), 
                        "pop_product" = population_weight[match(dist_counties %>% colnames(),
                                                         population_weight$CountyName), ]$weight * county_pop, 
                        "pop_diff" = population_weight[match(dist_counties %>% colnames(),
                                                             population_weight$CountyName), ]$weight - county_pop
                        ) 
  
  # assign county names and filter out only positive values to avoid double counting
  diff_all_counties <- rbind(diff_all_counties, 
                             left_join(diff_cases_df, 
                                       dist_df, 
                                       by = "CountyName")
                             ) %>% filter(cases >= 0)
}

# compute mean difference in 1-lag COVID-19 ID across time for each weighting 
# scheme
mean_diff_all_counties <- diff_all_counties %>% 
  group_by(CountyName, dist, pop_product, pop_diff) %>% 
  summarise(cases = mean(cases)) %>% 
  ungroup()  

# visualise the log transformed distance between counties against the 
# difference in 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = dist %>% log(), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log distance between counties") +
  ylab("Average difference in COVID-19 ID")

# visualise the log transformed INV-D weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(1 / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log inverse distance between counties") +
  ylab("Average difference in COVID-19 ID")

# visualise the log transformed PB weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(pop_product / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log population-based weights") +
  ylab("Average difference in COVID-19 ID")

# visualise the log transformed PD weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(abs(pop_diff) / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log abs. difference in population sizes times inverse distance") +
  ylab("Average difference in COVID-19 ID")

# Benchmark ARIMA  --------------------------------------------------------
# construct ARIMA models with optimal p, q for each county in isolation 
results_arima <- list()
for (county in covid_cases %>% colnames()) {
  
  covid_cases_county <- covid_cases_df %>% 
    dplyr::select(county %>% all_of()) %>% 
    ts(frequency = 52, start = c(2020, 1))
  
  arima_model <- auto.arima(y = covid_cases_county, 
                            d = 0)

  results_arima[[county]] <- list("model" = arima_model, 
                                  data.frame("BIC" = arima_model %>% BIC(), 
                                             "AIC" = arima_model %>% AIC()), 
                                  "arma" = paste0(arima_model$arma, 
                                                  collapse = "-"))
}

mean_results_arima <- lapply(results_arima, FUN = function(i) {
  i[[2]]
}) %>% 
  list.rbind() %>% 
  summarise(mean_BIC = mean(BIC), 
            mean_AIC = mean(AIC))
# p - q - d
lapply(results_arima, FUN = function(i) {
  i[[3]]
}) %>% list.rbind() %>% table()


# predict and compute MASE for ARIMA models 
mase_arima <- fit_and_predict_arima(forecast_window = 10)


# Queen's contiguity ------------------------------------------------------
# fit GNAR models for INV-D weighting 
results_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                          county_index = county_index_queen, 
                                          inverse_distance = TRUE)
return_best_model(results_queen)

# fit GNAR models for PB weighting 
results_pop_weighted_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                       county_index = county_index_queen,
                                                       inverse_distance = FALSE, 
                                                       weight_index = population_weight)
return_best_model(results_pop_weighted_queen)


# fit GNAR model for INV-D weighting with classification 
results_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                             county_index = county_index_queen,
                                                             inverse_distance = TRUE, 
                                                             weight_index = population_weight, 
                                                             weight_factor = urbanisation_factor,
                                                             globalalpha = TRUE)
return_best_model(results_class_queen)

# fit GNAR model for PB weighting with classification 
results_pop_weighted_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                             county_index = county_index_queen,
                                                             inverse_distance = FALSE, 
                                                             weight_index = population_weight, 
                                                             weight_factor = urbanisation_factor,
                                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_queen)


# fit GNAR models for SPL weighting 
results_old_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                              old = TRUE)
return_best_model(results_old_queen)

# fit GNAR models for SPL weighting with classification 
results_old_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                              old = TRUE, 
                                              weight_factor = urbanisation_factor,
                                              globalalpha = TRUE)
return_best_model(results_old_class_queen)

# Economic hub ------------------------------------------------------------
# fit GNAR models for INV-D weighting 
results_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                             county_index = county_index_eco_hubs, 
                                             inverse_distance = TRUE, 
                                             numeric_vertices = TRUE)
return_best_model(results_eco_hubs)

# fit GNAR models for PB weighting 
results_pop_weighted_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                          county_index = county_index_eco_hubs, 
                                                          inverse_distance = FALSE, 
                                                          numeric_vertices = TRUE)
return_best_model(results_pop_weighted_eco_hubs)

# fit GNAR models for INV-D weighting with classification 
results_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                   county_index = county_index_eco_hubs, 
                                                   inverse_distance = TRUE, 
                                                   numeric_vertices = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_eco_hubs)


# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                   county_index = county_index_eco_hubs, 
                                                   inverse_distance = FALSE, 
                                                   numeric_vertices = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_pop_weighted_class_eco_hubs)

# fit GNAR models for SPL weighting  
results_old_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                 old = TRUE)
return_best_model(results_old_eco_hubs)

# fit GNAR models for SPL weighting with classification 
results_old_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                 old = TRUE, 
                                                 weight_factor = urbanisation_factor, 
                                                 globalalpha = TRUE)
return_best_model(results_old_class_eco_hubs)

# Rail-based --------------------------------------------------------------
# fit GNAR models for INV-D weighting 
results_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                          county_index = county_index_train, 
                                          inverse_distance = TRUE, 
                                          numeric_vertices = TRUE)
return_best_model(results_train)

# fit GNAR models for PB weighting 
results_pop_weighted_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                          county_index = county_index_train, 
                                          inverse_distance = FALSE, 
                                          numeric_vertices = TRUE)
return_best_model(results_pop_weighted_train)

# fit GNAR models for INV-D weighting with classification  
results_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                county_index = county_index_train, 
                                                inverse_distance = TRUE, 
                                                numeric_vertices = TRUE, 
                                                weight_factor = urbanisation_factor, 
                                                globalalpha = TRUE)
return_best_model(results_class_train)

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                             county_index = county_index_train, 
                                                             inverse_distance = FALSE, 
                                                             numeric_vertices = TRUE, 
                                                             weight_factor = urbanisation_factor, 
                                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_train)

# fit GNAR models for SPL weighting 
results_old_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                              old = TRUE)
return_best_model(results_old_train)

# fit GNAR models for SPL weighting with classification 
results_old_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                    old = TRUE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
return_best_model(results_old_class_train)

# Delaunay triangulation --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                             county_index = county_index_delaunay, 
                                             inverse_distance = TRUE)
return_best_model(results_delaunay)

# fit GNAR models for PB weighting 
results_pop_weighted_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                          county_index = county_index_delaunay, 
                                                          inverse_distance = FALSE)
return_best_model(results_pop_weighted_delaunay)

# fit GNAR models for INV-D weighting with classification 
results_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                   county_index = county_index_delaunay, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_delaunay)

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                                county_index = county_index_delaunay, 
                                                                inverse_distance = FALSE, 
                                                                weight_factor = urbanisation_factor,
                                                                globalalpha = TRUE)
return_best_model(results_pop_weighted_class_delaunay)

# fit GNAR models for SPL weighting 
results_old_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                              old = TRUE)
return_best_model(results_old_delaunay)

# fit GNAR models for SPL weighting with classification 
results_old_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                       old = TRUE, 
                                                       weight_factor = urbanisation_factor, 
                                                       globalalpha = TRUE)
return_best_model(results_old_class_delaunay)

# Gabriel neighbourhood  --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = TRUE)
return_best_model(results_gabriel)

# fit GNAR models for PB weighting 
results_pop_weighted_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = FALSE)
return_best_model(results_pop_weighted_gabriel)

# fit GNAR models for INV-D weighting with classification 
results_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = TRUE, 
                                            weight_factor = urbanisation_factor, 
                                            globalalpha = TRUE)
return_best_model(results_class_gabriel)

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                  county_index = county_index_gabriel, 
                                                  inverse_distance = FALSE, 
                                                  weight_factor = urbanisation_factor, 
                                                  globalalpha = TRUE)
return_best_model(results_pop_weighted_class_gabriel)

# fit GNAR models for SPL weighting 
results_old_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                 old = TRUE)
return_best_model(results_old_gabriel)

# fit GNAR models for SPL weighting with classification 
results_old_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                      old = TRUE, 
                                                      weight_factor = urbanisation_factor, 
                                                      globalalpha = TRUE)
return_best_model(results_old_class_gabriel)

# Relative neighbourhood --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = TRUE)
return_best_model(results_relative)

# fit GNAR models for PB weighting 
results_pop_weighted_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = FALSE)
return_best_model(results_pop_weighted_relative)

# fit GNAR models for INV-D weighting with classification  
results_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = TRUE, 
                                             weight_factor = urbanisation_factor, 
                                             globalalpha = TRUE)
return_best_model(results_class_relative)

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = FALSE, 
                                             weight_factor = urbanisation_factor, 
                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_relative)

# fit GNAR models for SPL weighting 
results_old_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                                 old = TRUE)
return_best_model(results_old_relative)

# fit GNAR models for SPL weighting with classification 
results_old_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                                 old = TRUE, 
                                                 weight_factor = urbanisation_factor, 
                                                 globalalpha = TRUE)
return_best_model(results_old_class_relative)

# Sphere of influence neighbourhood ---------------------------------------
# fit GNAR models for INV_D weighting 
results_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = TRUE)
return_best_model(results_soi)

# fit GNAR models for PB weighting 
results_pop_weighted_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = FALSE)
return_best_model(results_pop_weighted_soi)

# fit GNAR models for INV-D weighting with classification  
results_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = TRUE, 
                                        weight_factor = urbanisation_factor, 
                                        globalalpha = TRUE)
return_best_model(results_class_soi)

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = FALSE, 
                                        weight_factor = urbanisation_factor, 
                                        globalalpha = TRUE)
return_best_model(results_pop_weighted_class_soi)

# fit GNAR models for SPL weighting 
results_old_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                                 old = TRUE)
return_best_model(results_old_soi)

# fit GNAR models for SPL weighting with classification 
results_old_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                            old = TRUE, 
                                            weight_factor = urbanisation_factor, 
                                            globalalpha = TRUE)
return_best_model(results_old_class_soi)


# Complete ----------------------------------------------------------------
results_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                             county_index = county_index_complete, 
                                             inverse_distance = TRUE)

return_best_model(results_complete)

# fit GNAR models for PB weighting
results_pop_weighted_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                          county_index = county_index_complete, 
                                                          inverse_distance = FALSE)

return_best_model(results_pop_weighted_complete)

# fit GNAR models for INV-D weighting with classification 
results_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                   county_index = county_index_complete, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_complete)

# fit GNAR models for PB weighting with classification  
results_pop_weighted_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                                county_index = county_index_complete, 
                                                                inverse_distance = FALSE, 
                                                                weight_factor = urbanisation_factor, 
                                                                globalalpha = TRUE)
return_best_model(results_pop_weighted_class_complete)

# fit GNAR models for SPL weighting 
results_old_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                 county_index = county_index_complete, 
                                                 old = TRUE)
return_best_model(results_old_complete)

# fit GNAR models for SPL weighting with classification 
results_old_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                       county_index = county_index_complete, 
                                                       old = TRUE, 
                                                       weight_factor = urbanisation_factor, 
                                                       globalalpha = TRUE)
return_best_model(results_old_class_complete)


# KNN ---------------------------------------------------------------------
# up to fully connected 
max_k <- dim(covid_cases)[2] - 1

# create list to save best performing model for each neighbourhood size 
knn_best <- list()

for (k in seq(1, max_k, by = 2)) {
  # create nb list for neighbourhood size k
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
  
  # create ordered coutny index data frame 
  county_index_knn <- data.frame("CountyName" = covid_net_knn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  # fit GNAR models for INV-D weighting 
  results_knn_id <- fit_and_predict_for_many(net = covid_net_knn, 
                                             county_index = county_index_knn, 
                                             inverse_distance = TRUE)
  # fit GNAR models for PB weighting 
  results_knn_pop <- fit_and_predict_for_many(net = covid_net_knn, 
                                              county_index = county_index_knn, 
                                              inverse_distance = FALSE)
  
  # fit GNAR models for INV-D weighting with classification
  results_knn_id_urban <- fit_and_predict_for_many(net = covid_net_knn, 
                                                   county_index = county_index_knn, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
  
  # fit GNAR models for PB weighting with classification 
  results_knn_pop_urban <- fit_and_predict_for_many(net = covid_net_knn, 
                                                    county_index = county_index_knn, 
                                                    inverse_distance = FALSE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # # fit GNAR models for SPL weighting 
  results_knn_old <- fit_and_predict_for_many(net = covid_net_knn, 
                                              county_index = county_index_knn, 
                                              old = TRUE)
  
  # fit GNAR models for SPL weighting with classification 
  results_knn_old_class <- fit_and_predict_for_many(net = covid_net_knn, 
                                                    county_index = county_index_knn, 
                                                    old = TRUE, 
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  
  # save best performing model for every k
  knn_best[[length(knn_best) + 1]] <- data.frame("k" = k, 
                                                 "best model id" = return_best_model(results_knn_id), 
                                                 "best model pop" = return_best_model(results_knn_pop),
                                                 "best model id class" = return_best_model(results_knn_id_urban), 
                                                 "best model pop class" = return_best_model(results_knn_pop_urban), 
                                                 "best model old" = return_best_model(results_knn_old), 
                                                 "best model old class" = return_best_model(results_knn_old_class)
  )
  
}

knn_best_df <- do.call(rbind.data.frame, knn_best)

# id - id.class - pop - pop.class - old - old.class
knn_best_res <- return_best_knn_dnn(knn_best_df)


knn_best_res$best.model.id.BIC <- knn_best_res$best.model.id.BIC %>% as.numeric()
knn_best_res$weighting <- c("INV-D", 
                            "INV-D+class",
                            "PB", 
                            "PB+class", 
                            "SPL", 
                            "SPL+class")


# DNN ---------------------------------------------------------------------
# find thresholds for DNN network 
max_dist_threshold <- dist_urbanisation %>% colMax() %>% max()
dist_urbanisation %>% colMax() %>% which.max()

min_dist_threshold <- dist_urbanisation %>% colMin() %>% max()

dist_urbanisation %>% colMin() %>% mean() # 42km
dist_urbanisation %>% colMin() %>% sd() # 18km

# create list to save best performing model for each distance threshold d
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
  
  # fit GNAR models for INV-D weighting 
  results_dnn_id <- fit_and_predict_for_many(net = covid_net_dnn, 
                                             county_index = county_index_dnn, 
                                             inverse_distance = TRUE)
  
  # fit GNAR models for PB weighting
  results_dnn_pop <- fit_and_predict_for_many(net = covid_net_dnn, 
                                              county_index = county_index_dnn, 
                                              inverse_distance = FALSE)
  
  # fit GNAR models for INV-D weighting with classification 
  results_dnn_id_urban <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                   county_index = county_index_dnn, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
  
  # fit GNAR models for PB weighting with classification 
  results_dnn_pop_urban <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                    county_index = county_index_dnn, 
                                                    inverse_distance = FALSE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # fit GNAR models for SPL weighting
  results_dnn_old <- fit_and_predict_for_many(net = covid_net_dnn, 
                                              county_index = county_index_dnn, 
                                              old = TRUE)
  
  # fit GNAR models for SPL weighting with classification 
  results_dnn_old_class <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                    county_index = county_index_dnn, 
                                                    old = TRUE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # save best performing model for every d
  dnn_best[[length(dnn_best) + 1]] <- data.frame("d in km" = d, 
                                                 "best model id" = return_best_model(results_dnn_id), 
                                                 "best model pop" = return_best_model(results_dnn_pop),
                                                 "best model id class" = return_best_model(results_dnn_id_urban), 
                                                 "best model pop class" = return_best_model(results_dnn_pop_urban),
                                                 "best model old" = return_best_model(results_dnn_old), 
                                                 "best model old class" = return_best_model(results_dnn_old_class)
  )
}

dnn_best_df <- do.call(rbind.data.frame, dnn_best)

# id - id.class - pop - pop.class - old - old.class 
dnn_best_res <- return_best_knn_dnn(dnn_best_df)
dnn_best_res$best.model.id.BIC <- dnn_best_res$best.model.id.BIC %>% as.numeric()

dnn_best_res$weighting <- c("INV-d", 
                            "INV-D+class",
                            "PB", 
                            "PB+class", 
                            "SPL", 
                            "SPL+class")




# Compare BIC for GNAR models ---------------------------------------------
# for latex 
bic_development <- cbind(results_delaunay[, c(1, 3)], 
                         results_gabriel[, 3], 
                         results_soi[, 3], 
                         results_relative[, 3], 
                         results_queen[, 3], 
                         results_eco_hubs[, 3], 
                         results_train[, 3])


strCaption <- "\\code{GNAR} models with INV-D weighting and without vertex classification
for \\textbf{Delaunay triangulation}, \\textbf{Gabriel}, \\textbf{SOI}, 
\\textbf{Relative neighbourhood}, \\textbf{Queen's contiguity}, 
\\textbf{Economic hub} and \\textbf{Railway-based} network and the corresponding
BIC for model selection; the best performing model for each network is 
highlighted in bold red"
print(xtable(bic_development,
             digits=2,
             caption=strCaption,
             label="tab:gnar_all", 
             align = c("", "l", "|", "r", "r", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(bic_development)),
                        command = c(paste("\\toprule \n",
                                          " Model name & Delaunay & Gabriel & 
                                          SOI & Relative & Queen & Eco. hub & 
                                          Railway \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# Performance overview ----------------------------------------------------
compare_df <- data.frame("network" = c(rep("Delaunay", 6),
                                       rep("Gabriel", 6),
                                       rep("SOI", 6),
                                       rep("Relative", 6), 
                                       rep("Complete", 6), 
                                       rep("Queen", 6),
                                       rep("Eco. hub", 6), 
                                       rep("Rail", 6)
                                       ), 
                         "hyperparameters" = rep(c("INV-D", 
                                                   "INV-D+class", 
                                                   "PB", 
                                                   "PB+class", 
                                                   "SPL", 
                                                   "SPL+class"), 
                                                 8), 
                         "model" = c(
                                     return_best_model(results_delaunay)$name, 
                                     return_best_model(results_class_delaunay)$name, 
                                     return_best_model(results_pop_weighted_delaunay)$name,
                                     return_best_model(results_pop_weighted_class_delaunay)$name,
                                     return_best_model(results_old_delaunay)$name,
                                     return_best_model(results_old_class_delaunay)$name,
                                     
                                     return_best_model(results_gabriel)$name, 
                                     return_best_model(results_class_gabriel)$name, 
                                     return_best_model(results_pop_weighted_gabriel)$name,
                                     return_best_model(results_pop_weighted_class_gabriel)$name,
                                     return_best_model(results_old_gabriel)$name,
                                     return_best_model(results_old_class_gabriel)$name,
                                     
                                     return_best_model(results_soi)$name, 
                                     return_best_model(results_class_soi)$name, 
                                     return_best_model(results_pop_weighted_soi)$name,
                                     return_best_model(results_pop_weighted_class_soi)$name, 
                                     return_best_model(results_old_soi)$name,
                                     return_best_model(results_old_class_soi)$name,
                                     
                                     return_best_model(results_relative)$name, 
                                     return_best_model(results_class_relative)$name, 
                                     return_best_model(results_pop_weighted_relative)$name,
                                     return_best_model(results_pop_weighted_class_relative)$name,
                                     return_best_model(results_old_relative)$name,
                                     return_best_model(results_old_class_relative)$name,
                                     
                                     return_best_model(results_complete)$name, 
                                     return_best_model(results_class_complete)$name, 
                                     return_best_model(results_pop_weighted_complete)$name,
                                     return_best_model(results_pop_weighted_class_complete)$name,
                                     return_best_model(results_old_complete)$name, 
                                     return_best_model(results_old_class_complete)$name, 
                                     
                                     return_best_model(results_queen)$name, 
                                     return_best_model(results_class_queen)$name, 
                                     return_best_model(results_pop_weighted_queen)$name,
                                     return_best_model(results_pop_weighted_class_queen)$name,
                                     return_best_model(results_old_queen)$name,
                                     return_best_model(results_old_class_queen)$name, 
                                     
                                     return_best_model(results_eco_hubs)$name, 
                                     return_best_model(results_class_eco_hubs)$name, 
                                     return_best_model(results_pop_weighted_eco_hubs)$name,
                                     return_best_model(results_pop_weighted_class_eco_hubs)$name,
                                     return_best_model(results_old_eco_hubs)$name,
                                     return_best_model(results_old_class_eco_hubs)$name,
                                     
                                     return_best_model(results_train)$name, 
                                     return_best_model(results_class_train)$name, 
                                     return_best_model(results_pop_weighted_train)$name,
                                     return_best_model(results_pop_weighted_class_train)$name,
                                     return_best_model(results_old_train)$name,
                                     return_best_model(results_old_class_train)$name
                                     ), 
                         "BIC" = c(
                                   return_best_model(results_delaunay)$BIC, 
                                   return_best_model(results_class_delaunay)$BIC, 
                                   return_best_model(results_pop_weighted_delaunay)$BIC,
                                   return_best_model(results_pop_weighted_class_delaunay)$BIC,
                                   return_best_model(results_old_delaunay)$BIC,
                                   return_best_model(results_old_class_delaunay)$BIC,
                                   
                                   return_best_model(results_gabriel)$BIC, 
                                   return_best_model(results_class_gabriel)$BIC, 
                                   return_best_model(results_pop_weighted_gabriel)$BIC,
                                   return_best_model(results_pop_weighted_class_gabriel)$BIC,
                                   return_best_model(results_old_gabriel)$BIC,
                                   return_best_model(results_old_class_gabriel)$BIC,
                                   
                                   return_best_model(results_soi)$BIC, 
                                   return_best_model(results_class_soi)$BIC, 
                                   return_best_model(results_pop_weighted_soi)$BIC,
                                   return_best_model(results_pop_weighted_class_soi)$BIC, 
                                   return_best_model(results_old_soi)$BIC,
                                   return_best_model(results_old_class_soi)$BIC,
                                   
                                   return_best_model(results_relative)$BIC, 
                                   return_best_model(results_class_relative)$BIC, 
                                   return_best_model(results_pop_weighted_relative)$BIC,
                                   return_best_model(results_pop_weighted_class_relative)$BIC,
                                   return_best_model(results_old_relative)$BIC,
                                   return_best_model(results_old_class_relative)$BIC,
                                   
                                   return_best_model(results_complete)$BIC, 
                                   return_best_model(results_class_complete)$BIC, 
                                   return_best_model(results_pop_weighted_complete)$BIC,
                                   return_best_model(results_pop_weighted_class_complete)$BIC, 
                                   return_best_model(results_old_complete)$BIC, 
                                   return_best_model(results_old_class_complete)$BIC,
                                   
                                   return_best_model(results_queen)$BIC, 
                                   return_best_model(results_class_queen)$BIC, 
                                   return_best_model(results_pop_weighted_queen)$BIC,
                                   return_best_model(results_pop_weighted_class_queen)$BIC,
                                   return_best_model(results_old_queen)$BIC, 
                                   return_best_model(results_old_class_queen)$BIC,
                                   
                                   return_best_model(results_eco_hubs)$BIC, 
                                   return_best_model(results_class_eco_hubs)$BIC, 
                                   return_best_model(results_pop_weighted_eco_hubs)$BIC,
                                   return_best_model(results_pop_weighted_class_eco_hubs)$BIC,
                                   return_best_model(results_old_eco_hubs)$BIC, 
                                   return_best_model(results_old_class_eco_hubs)$BIC, 
                                   
                                   return_best_model(results_train)$BIC, 
                                   return_best_model(results_class_train)$BIC, 
                                   return_best_model(results_pop_weighted_train)$BIC,
                                   return_best_model(results_pop_weighted_class_train)$BIC,
                                   return_best_model(results_old_train)$BIC,
                                   return_best_model(results_old_class_train)$BIC
                         ))

# for latex 
strCaption <- "Comparison of best performing \\code{GNARfit()} models for 
COVID-19 networks, for all weighting schemes with and without classification 
(+class) according to urbanisation of county; 
the best performing model is highlighted in bold red"
print(xtable(compare_df,
             digits=2,
             caption=strCaption,
             label="tab:overview_weighting_best_model", 
             align = c("", "l", "|", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(compare_df)),
                        command = c(paste("\\toprule \n",
                                          " Network & hyperparameter & \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# for KNN and DNN
compare_knn_dnn <- data.frame("network" = c(rep("KNN", 6), 
                                            rep("DNN", 6)), 
                              "hyperparameters" = rep(c("INV-D", 
                                                        "INV-D+class", 
                                                        "PB", 
                                                        "PB+class", 
                                                        "SPL", 
                                                        "SPL+class"), 
                                                      2), 
                              "k_d" = c(knn_best_res$k, 
                                        dnn_best_res$d.in.km), 
                              "model" = c(knn_best_res$best.model.id.name, 
                                          dnn_best_res$best.model.id.name), 
                              "BIC" = c(knn_best_res$best.model.id.BIC, 
                                        dnn_best_res$best.model.id.BIC) %>% 
                                as.numeric() %>% 
                                round(2)
                              )


# for latex 
strCaption <- "Comparison of best performing \\code{GNARfit()} models for KNN 
and DNN networks, for all weighting schemes with and without classification 
(+class) according to urbanisation of county; 
the best performing model is highlighted in bold red"
print(xtable(compare_knn_dnn,
             digits=2,
             caption=strCaption,
             label="tab:overview_best_model_knn_dnn", 
             align = c("", "l", "|", "r", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(compare_knn_dnn)),
                        command = c(paste("\\toprule \n",
                                          " Network & hyperparameter & k / d [km] & 
                                          \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# best KNN / DNN 
compare_knn_dnn %>% 
  group_by(network) %>% 
  summarise(which.min(BIC))
# first model for both KNN / DNN is best performing 

# absolute best for KNN / DNN 
compare_knn_dnn[which.min(compare_knn_dnn$BIC), ]

# absolute best 
compare_df[which.min(compare_df$BIC), ]


# Best KNN ----------------------------------------------------------------
# construct best performing KNN network and GNAR model
# create nb list
opt_knn_net <- knearneigh(x = coord_urbanisation, 
                          k = 21, 
                          longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames())

# create igraph object
opt_knn_net_igraph <- neighborsDataFrame(nb = opt_knn_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
opt_knn_net_gnar <- opt_knn_net_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame
county_index_opt_knn <- data.frame("CountyName" = opt_knn_net_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# compute network characteristics 
graph_char_knn <- network_characteristics(opt_knn_net_igraph, 
                                          "KNN")


# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(opt_knn_net,
#      coord_urbanisation,
#      add = TRUE,
#      pch = 19, cex = 0.6)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

opt_knn_mod <- fit_and_predict(net = opt_knn_net_gnar, 
                               alpha = 5, 
                               beta = c(1, 1, 1, 1, 0), 
                               globalalpha = TRUE, 
                               old = TRUE, 
                               county_index = county_index_opt_knn, 
                               return_model = TRUE)
summary(opt_knn_mod$mod)
BIC(opt_knn_mod)
AIC(opt_knn_mod)
opt_knn_mod %>% coef()


# compute Moran's I
moran_knn <- moran_I_permutation_test(data = COVID_weekly_data, 
                                      g = opt_knn_net_igraph, 
                                      name = "knn")
# # visualise Moran's I
# ggplot(moran_knn, 
#        aes(x = dates, 
#            y = moran)) +
#   geom_line() +
#   xlab("Time") +
#   ylab("Moran's I") +
#   geom_vline(aes(xintercept = as.Date("18.08.2020",
#                                       format = "%d.%m.%Y"), 
#                  color = "County-specific restrictions")) +
#   geom_vline(aes(xintercept = as.Date("26.12.2020", 
#                                       format = "%d.%m.%Y"), 
#                  color = "Level-5 lockdown")) +
#   geom_vline(aes(xintercept = as.Date("26.07.2021",
#                                       format = "%d.%m.%Y"), 
#                  color = "Indoor dining")) +
#   geom_vline(aes(xintercept = as.Date("06.03.2022",
#                                       format = "%d.%m.%Y"), 
#                  color = "End")) +
#   scale_color_brewer(palette = "Set1") + 
#   theme(legend.position = "None")
# ggsave("Figures/MoransI/covid_moran_knn.pdf", 
#        width = 27, height = 14, unit = "cm")


# Best DNN ----------------------------------------------------------------
# construct best performing DNN network and GNAR model 
# create nb list
opt_dnn_net <- dnearneigh(x = coord_urbanisation, 
                          d1 = 0, 
                          d2 = 325,
                          row.names = coord_urbanisation %>% rownames(),
                          longlat = TRUE, 
                          use_s2 = TRUE) 

# create igraph object 
opt_dnn_net_igraph <- neighborsDataFrame(opt_dnn_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

# create GNAR object 
opt_dnn_net_gnar <- opt_dnn_net_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame
county_index_opt_dnn <- data.frame("CountyName" = opt_dnn_net_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# compute network characteristics
graph_char_dnn <- network_characteristics(opt_dnn_net_igraph, 
                                          "DNN")



# for latex 
strCaption <- paste0("Summary for the \\textbf{DNN} network, av. short for 
                     average, s.d. short for standard deviation")
print(xtable(graph_char_dnn,
             digits=2,
             caption=strCaption,
             label="tab:summary_dnn", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_char_dnn)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# for complete graph
graph_complete <- network_characteristics(complete_net_igraph, 
                                          "complete")

# Network characteristics for KNN and DNN
graph_knn_dnn <- cbind(graph_char_knn, 
                       graph_char_dnn$DNN, 
                       graph_complete$complete)

strCaption <- "Overview of network characteristics for best performing KNN (k = 21),
DNN (d = 175) and Complete network, including average (av.) degree, density, average (av.) 
shortest path length (SPL), global and average (av.) local clustering (clust.) 
as well as average (av.) betweenness (betw.) and its standard deviation (s.d.)"
print(xtable(graph_knn_dnn,
             digits=2,
             caption=strCaption,
             label="tab:summary_knn_dnn", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_knn_dnn)),
                        command = c(paste("\\toprule \n",
                                          " Metric & KNN & DNN & Complete \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(opt_dnn_net,
#      coord_urbanisation,
#      add = TRUE,
#      pch = 19, cex = 0.6)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)


# best model 
opt_dnn_mod <- fit_and_predict(net = opt_dnn_net_gnar, 
                               alpha = 5, 
                               beta = c(1, 1, 1, 1, 0), 
                               globalalpha = TRUE, 
                               inverse_distance = TRUE, 
                               county_index = county_index_opt_dnn, 
                               return_model = TRUE)
summary(opt_dnn_mod$mod)
BIC(opt_dnn_mod)


# compute Moran's I
moran_dnn <- moran_I_permutation_test(data = COVID_weekly_data, 
                                      g = opt_dnn_net_igraph, 
                                      name = "dnn")

# visualise Moran's I
# ggplot(moran_dnn, 
#        aes(x = dates, 
#            y = moran)) +
#   geom_line() +
#   xlab("Time") +
#   ylab("Moran's I") +
#   geom_vline(aes(xintercept = as.Date("18.08.2020",
#                                       format = "%d.%m.%Y"), 
#                  color = "County-specific restrictions")) +
#   geom_vline(aes(xintercept = as.Date("26.12.2020", 
#                                       format = "%d.%m.%Y"), 
#                  color = "Level-5 lockdown")) +
#   geom_vline(aes(xintercept = as.Date("26.07.2021",
#                                       format = "%d.%m.%Y"), 
#                  color = "Indoor dining")) +
#   geom_vline(aes(xintercept = as.Date("06.03.2022",
#                                       format = "%d.%m.%Y"), 
#                  color = "End")) +
#   scale_color_brewer(palette = "Set1") + 
#   theme(legend.position = "None")
# ggsave("Figures/MoransI/covid_moran_dnn.pdf", 
#        width = 27, height = 14, unit = "cm")

# Best model for each network ---------------------------------------------
# find best model for each sparse and Complete network 
best_model <- compare_df %>% 
  group_by(network) %>% 
  summarise(BIC = min(BIC))

# filter minimum BIC for KNN and DNN network 
best_model_knn_dnn <- compare_knn_dnn %>% 
  group_by(network) %>% 
  summarise(BIC = min(BIC))


best_model_overview <- left_join(best_model, 
                                 compare_df, 
                                 by = c("BIC", "network")) %>% 
  dplyr::select(network, 
                hyperparameters, 
                model, 
                BIC) %>% 
  as.data.frame()

best_model_knn_dnn_overview <- left_join(best_model_knn_dnn, 
                                         compare_knn_dnn, 
                                         by = c("BIC", "network")) %>% 
  dplyr::select(network, 
                hyperparameters, 
                model, 
                BIC) %>% 
  as.data.frame()


# merge data frames with best models for sparse + Complete networks and 
# KNN + DNN networks 
best_model_overview_all <- rbind(best_model_overview, 
                                 best_model_knn_dnn_overview)


# Residual analysis  ------------------------------------------------------

# for each network, fit best performing model and analysis residuals 
# in scatter plot and qqplot 

# Queen 
model_queen <- fit_and_predict(alpha = 4, 
                               beta = c(2, 1, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_queen_gnar, 
                               old = TRUE, 
                               county_index = county_index_queen, 
                               return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_queen <- check_and_plot_residuals(model = model_queen, 
                                            network_name = "queen")


# Eco hubs
model_eco_hubs <- fit_and_predict(alpha = 5, 
                                  beta = c(2, 2, 1, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_eco_hubs_gnar, 
                                  old = TRUE, 
                                  numeric_vertices = TRUE, 
                                  county_index = county_index_eco_hubs,
                                  return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_eco_hubs <- check_and_plot_residuals(model = model_eco_hubs, 
                                               network_name = "eco_hubs", 
                                               alpha = 5)


# Train 
model_train <- fit_and_predict(alpha = 4, 
                               beta = c(2, 1, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_train_gnar, 
                               old = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train,
                               return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_train <- check_and_plot_residuals(model = model_train, 
                                            network_name = "train")


# Delaunay 
model_delaunay <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 2, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_delaunay_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_delaunay, 
                                  return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_delaunay <- check_and_plot_residuals(model = model_delaunay, 
                                               network_name = "delaunay")

# Gabriel
model_gabriel <- fit_and_predict(alpha = 4, 
                                 beta = c(2, 1, 1, 1),
                                 globalalpha = TRUE, 
                                 net = covid_net_gabriel_gnar, 
                                 old = TRUE, 
                                 county_index = county_index_gabriel,
                                 return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_gabriel <- check_and_plot_residuals(model = model_gabriel, 
                                              network_name = "gabriel", 
                                              alpha = 4)


# Relative 
model_relative <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 1, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_relative_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_relative,
                                  return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_relative <- check_and_plot_residuals(model = model_relative, 
                                               network_name = "relative")

# SOI
model_soi <- fit_and_predict(alpha = 4, 
                             beta = c(2, 1, 1, 1),
                             globalalpha = TRUE, 
                             net = covid_net_soi_gnar, 
                             old = TRUE, 
                             county_index = county_index_soi, 
                             return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_soi <- check_and_plot_residuals(model = model_soi, 
                                         network_name = "soi")

# KNN 
# compute residuals for fitted values and plot in scatter plot
residuals_knn <- check_and_plot_residuals(model = opt_knn_mod, 
                                          network_name = "KNN", 
                                          alpha = 5)

# DNN
# compute residuals for fitted values and plot in scatter plot
residuals_dnn <- check_and_plot_residuals(model = opt_dnn_mod, 
                                          network_name = "DNN", 
                                          alpha = 5)

# Complete 
model_complete <- fit_and_predict(alpha = 5,
                                  beta = c(1, 1, 1, 1, 0), 
                                  globalalpha = TRUE, 
                                  net = complete_net_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_complete,
                                  return_model = TRUE)

# compute residuals for fitted values and plot in scatter plot
residuals_complete <- check_and_plot_residuals(model = model_complete, 
                                               network_name = "complete", 
                                               alpha = 5)


# add AIC 
best_model_overview_all$AIC <- c(AIC(model_complete),
                                 AIC(model_delaunay), 
                                 AIC(model_eco_hubs), 
                                 AIC(model_gabriel), 
                                 AIC(model_queen), 
                                 AIC(model_train), 
                                 AIC(model_relative), 
                                 AIC(model_soi), 
                                 AIC(opt_dnn_mod), 
                                 AIC(opt_knn_mod)
                                 )

# for latex 
strCaption <- "Overview over best performing \\code{GNARfit()} models for 
each COVID-19 network, SPL weighting best performing for all networks"
print(xtable(best_model_overview_all[c(2, 4, 8, 7, 1, 10, 9, 5, 3, 6), ],
             digits=2,
             caption=strCaption,
             label="tab:best_models", 
             align = c("", "l", "|", "r", "c", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_model_overview_all)),
                        command = c(paste("\\toprule \n",
                                          " Network & weighing & 
                                          \\code{GNAR} model & BIC & 
                                          AIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# MASE --------------------------------------------------------------------
# compute MASE for the last 10 weeks 

# Queen
mase_queen <- compute_MASE(model = model_queen, 
                           network_name = "Queen",
                           counties = counties)

# Eco hubs
mase_eco_hubs <- compute_MASE(model = model_eco_hubs, 
                              network_name = "Eco hubs",
                              counties = counties)

# Train 
mase_train <- compute_MASE(model = model_train, 
                           network_name = "Train",
                           counties = counties)

# Delaunay 
mase_delaunay <- compute_MASE(model = model_delaunay, 
                              network_name = "Delaunay",
                              counties = counties)

# Gabriel 
mase_gabriel <- compute_MASE(model = model_gabriel, 
                             network_name = "Gabriel", 
                             counties = counties)

# Relative 
mase_relative <- compute_MASE(model = model_relative, 
                              network_name = "Relative",
                              counties = counties)

# Soi
mase_soi <- compute_MASE(model = model_soi, 
                         network_name = "SOI",
                         counties = counties)

# KNN
mase_knn <- compute_MASE(model = opt_knn_mod, 
                         network_name = "KNN",
                         counties = counties)

# DNN
mase_dnn <- compute_MASE(model = opt_dnn_mod, 
                         network_name = "DNN",
                         counties = counties)

# Complete 
mase_complete <- compute_MASE(model = model_complete, 
                              network_name = "Complete",
                              counties = counties)

# summarise all MASE values for each network in a data frame  
mase_overview <- rbind.data.frame(mase_queen, 
                                  mase_eco_hubs, 
                                  mase_train, 
                                  mase_delaunay, 
                                  mase_gabriel, 
                                  mase_relative, 
                                  mase_soi, 
                                  mase_knn, 
                                  mase_dnn, 
                                  mase_complete, 
                                  mase_arima) %>% 
  na.omit()




# plot MASE for counties: Dublin, Wicklow, Kerry, Donegal
ggplot(mase_overview, 
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
  scale_color_manual(values = c("ARIMA" = "#3A3B3C", 
                                "KNN" = "#F8766D", 
                                "DNN" = "#D89000", 
                                "Complete" = "#A3A500", 
                                "Queen" = "#E76BF3", 
                                "Eco hubs" = "#00BF7D", 
                                "Gabriel" = "#00BFC4", 
                                "Relative" = "#00B0F6", 
                                "SOI" = "#9590FF", 
                                "Delaunay" = "#39B600", 
                                "Train" = "#FF62BC"), 
                     name = "Network")
ggsave("Figures/GNAR_entire_dataset/mase.pdf", 
       width = 26, height = 13, units = "cm")

# plot MASE for the Delaunay, Gabriel, Relative and SOI as well as Rail-way 
# based network 
ggplot(mase_overview %>% filter(type %in% c("Delaunay", 
                                            "Gabriel", 
                                            "Relative", 
                                            "SOI", 
                                            "Train", 
                                            "ARIMA")), 
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
  scale_color_manual(values = c("ARIMA" = "#3A3B3C", 
                                "Gabriel" = "#00BFC4", 
                                "Relative" = "#00B0F6", 
                                "SOI" = "#9590FF", 
                                "Delaunay" = "#39B600", 
                                "Train" = "#FF62BC"), 
                     name = "Network")
ggsave("Figures/GNAR_entire_dataset/mase_delaunay_etc_zoom.pdf", 
       width = 26, height = 13, units = "cm")

# plot MASE for the KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_overview %>% filter(type %in% c("DNN", 
                                            "KNN", 
                                            "Queen", 
                                            "Eco hubs", 
                                            "Complete", 
                                            "ARIMA")), 
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
  scale_color_manual(values = c("ARIMA" = "#3A3B3C", 
                                "KNN" = "#F8766D", 
                                "DNN" = "#D89000", 
                                "Complete" = "#A3A500", 
                                "Queen" = "#E76BF3", 
                                "Eco hubs" = "#00BF7D"), 
                     name = "Network")
ggsave("Figures/GNAR_entire_dataset/mase_knn_etc_zoom.pdf", 
       width = 26, height = 13, units = "cm")


# KS-test -----------------------------------------------------------------
ks_queen <- ks_residuals(mase_queen)
ks_eco <- ks_residuals(mase_eco_hubs)
ks_train <- ks_residuals(mase_train)
ks_delaunay <- ks_residuals(mase_delaunay)
ks_gabriel <- ks_residuals(mase_gabriel)
ks_relative <- ks_residuals(mase_relative)
ks_soi <- ks_residuals(mase_soi)
ks_complete <- ks_residuals(mase_complete)
ks_knn <- ks_residuals(mase_knn)
ks_dnn <- ks_residuals(mase_dnn)



# Scale-free --------------------------------------------------------------
# analyse log-log behaviour and regression R squared
queen_scale_free <- is_scale_free(igraph_net = covid_net_queen_igraph, 
                                  network_name = "queen")
queen_scale_free$graph
queen_scale_free$R_squared

eco_hub_scale_free <- is_scale_free(covid_net_eco_hubs_igraph, 
                                    network_name = "eco_hubs")
eco_hub_scale_free$graph
eco_hub_scale_free$R_squared

train_scale_free <- is_scale_free(covid_net_train_igraph, 
                                  network_name = "train")
train_scale_free$graph
train_scale_free$R_squared

knn_scale_free <- is_scale_free(opt_knn_net_igraph,
                                network_name = "knn")
knn_scale_free$graph
knn_scale_free$R_squared

dnn_scale_free <- is_scale_free(opt_dnn_net_igraph, 
                                network_name = "dnn")
dnn_scale_free$graph
dnn_scale_free$R_squared

delaunay_scale_free <- is_scale_free(covid_net_delaunay_igraph, 
                                     network_name = "delaunay")
delaunay_scale_free$graph
delaunay_scale_free$R_squared

gabriel_scale_free <- is_scale_free(covid_net_gabriel_igraph, 
                                    network_name = "gabriel")
gabriel_scale_free$graph
gabriel_scale_free$R_squared

relative_scale_free <- is_scale_free(covid_net_relative_igraph, 
                                     network_name = "relative")
relative_scale_free$graph
relative_scale_free$R_squared

soi_scale_free <- is_scale_free(igraph_net = covid_net_soi_igraph, 
                                network_name = "soi")
soi_scale_free$graph
soi_scale_free$R_squared

# for latex
scale_free_overview <- data.frame("name" = c("Delaunay triangulation", 
                                             "Gabriel", 
                                             "SOI",                                           
                                             "Relative neighbourhood", 
                                             "KNN", 
                                             "DNN",
                                             "Queen's", 
                                             "Economic hub", 
                                             "Railway-based"), 
                                  "R_squared" = c(delaunay_scale_free$R_squared, 
                                                  gabriel_scale_free$R_squared,
                                                  soi_scale_free$R_squared, 
                                                  relative_scale_free$R_squared,
                                                  knn_scale_free$R_squared,
                                                  dnn_scale_free$R_squared, 
                                                  queen_scale_free$R_squared, 
                                                  eco_hub_scale_free$R_squared, 
                                                  train_scale_free$R_squared))

strCaption <- "Test for scale-free property for each constructed network, 
$R^2$ for regression of log empirical cumulative distribution on log transformed 
degree"
print(xtable(scale_free_overview,
             digits=2,
             caption=strCaption,
             label="tab:scale_free", 
             align = c("", "l", "|", "c")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(scale_free_overview)),
                        command = c(paste("\\toprule \n",
                                          "Network & $R^2$ \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# Save KNN / DNN network --------------------------------------------------
save(list = c("opt_knn_net_igraph", 
              "opt_dnn_net_igraph"), 
     file = "Data/RObjects/KNN_DNN_igraph.RData")

save(list = c("opt_knn_net_gnar", 
              "opt_dnn_net_gnar"), 
     file = "Data/RObjects/KNN_DNN_GNAR.RData")




