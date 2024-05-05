### FUNCTIONS

# source adapted functions from GNAR package 
source("GNAR/GNARdesign.R")
source("GNAR/GNARfit.R")
source("GNAR/NofNeighbours.R")


squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

# Dataframe functions -----------------------------------------------------
# find maximum value for each column
colMax <- function(data) {
  return(sapply(data, max, na.rm = TRUE))
}

# find minimum value for each column 
colMin <- function(data) {
  return(sapply(data, min, na.rm = TRUE))
}

# Economic hub network ----------------------------------------------------
# create economic hub network based on Queen's contiguity network and edges to
# nearest economic hub 
create_eco_hub_igraph <- function(dist_df, # pairwise distance between vertices in data frame
                                  coord) { # matrix of vertex  coordinates
  queen_adj <- covid_net_queen_igraph %>% 
    as_adjacency_matrix()
  
  for (county in queen_adj %>% rownames()) {
    # select nearest economic hub for county in question 
    nearest_hub <- dist_df[county, hubs] %>% which.min() %>% names()
    
    queen_adj[county, nearest_hub] <- 1
    queen_adj[nearest_hub, county] <- 1
  }
  
  index_list <- data.frame("county" = queen_adj %>% rownames(), 
                           "index" = seq(1, 26))
  
  # assign numbers for row names and column names  
  rownames(queen_adj) <- seq(1, 26)
  colnames(queen_adj) <- seq(1, 26)
  
  # create igraph object from adjacency matrix
  covid_net_eco_hubs_igraph <- graph_from_adjacency_matrix(adjmatrix = queen_adj, 
                                                           mode = "undirected") 
  
  # order coordinates in the same order as adjacency matrix
  coord_ord <- coord[match(index_list$county, 
                           rownames(coord)), ] %>% 
    as.data.frame()
  coord_hubs <- coord_ord %>% filter(rownames(coord_ord) %in% hubs)
  
  return(list("igraph_net" = covid_net_eco_hubs_igraph, 
              "ordered_coord" = coord_ord, 
              "ordered_hubs" = coord_hubs))
}


# igraph ------------------------------------------------------------------
# generate nb list from igraph object
# requires numeric vertices (i.e. no county names for vertex names)
igraph2nb <- function(gr) {
  edges <- get.edgelist(gr)
  mode(edges) <- "integer"
  return(neig(edges = edges) %>% neig2nb())
}


neighborsDataFrame <- function(nb) {
  
  ks = data.frame(k = unlist(mapply(rep, 1:length(nb), 
                                    sapply(nb, length), 
                                    SIMPLIFY = FALSE) ), 
                  k_nb = unlist(nb) )
  
  nams = data.frame(id = attributes(nb)$region.id, 
                    k = 1:length(nb))
  
  o = merge(ks, nams, 
            by.x = 'k', 
            by.y = 'k')
  o = merge(o, nams, 
            by.x = 'k_nb', 
            by.y = 'k', 
            suffixes = c("","_neigh"))
  
  o[, c("id", "id_neigh")] %>% return()
}

# Network characteristics -------------------------------------------------
# compute network characteristics
network_characteristics <- function(igraph_obj, 
                                    network_name) {
  
  density <- igraph_obj %>% graph.density() 
  apl <- igraph_obj %>%  average.path.length(directed = FALSE) 
  
  global_clust <- igraph_obj %>% transitivity(type = "global") 
  mean_local_clust <- igraph_obj %>% 
    transitivity(type = "local",
                 isolates = "zero") %>% 
    mean()
  
  
  degree_v <- igraph_obj %>% igraph::degree()
  av_degree <- degree_v %>% mean() 
  
  # model Bernoulli Random Graph to check for small world behaviour 
  brg <- erdos.renyi.game(n = igraph_obj %>% gorder(), 
                          p.or.m = igraph_obj %>% gsize(), 
                          directed = FALSE,
                          loops = FALSE,
                          type = "gnm")
  
  apl_brg <- brg %>% average.path.length(directed = FALSE)
  mean_clustering_brg <- brg %>% 
    transitivity(type = "local",
                 isolates = "zero") %>% 
    mean()
  
  
  max_degree <- degree_v %>% max() 
  which(degree_v == max_degree) 
  
  min_degree <- degree_v %>% min() 
  which(degree_v == min_degree) 
  
  # betweenness
  bet <- betweenness(igraph_obj, 
                     v=V(igraph_obj), 
                     directed = FALSE)
  
  min_bet <- bet %>% min() 
  bet[which(bet == min_bet)] 
  
  max_bet <- bet %>% max()
  bet[which(bet == max_bet)] 
  
  
  graph_char <- data.frame("metric" = c("av. degree", 
                                        "density", 
                                        "av. SPL", 
                                        "global clust.", 
                                        "av. local clust.", 
                                        "av. betw.", 
                                        "s.d. betw.", 
                                        "BRG av. SPL", 
                                        "BRG av. local clust."), 
                           "values" = c(av_degree, 
                                        density, 
                                        apl, 
                                        global_clust, 
                                        mean_local_clust, 
                                        mean(bet), 
                                        sd(bet),
                                        apl_brg, 
                                        mean_clustering_brg)) 
  colnames(graph_char) <- c("metric", network_name)
  return(graph_char)
}

# Correlation -------------------------------------------------------------
# compute Residual Sum of Squares
RSS <- function(GNARfit_object) {
  return(sum(residuals(GNARfit_object)^2))
}

# Moran's I permutation test 
moran_I_permutation_test <- function(data = COVID_weekly_data,
                                     g, 
                                     county_index = NULL, 
                                     name, 
                                     time_col = "yw", 
                                     cases_col = "weeklyCases") {
  
  # compute shortest path length for each vertex pair
  distMatrix <- exp(shortest.paths(g, v=V(g), to=V(g)) * (-1))
  
  # if igraph does not have county names for vertices 
  if (!is.null(county_index)) {
    # assign names 
    county_ordering <- match(seq(1, 26), 
                             county_index$index)
    county_names <- county_index[county_ordering, ]$CountyName
    
    rownames(distMatrix) <- county_names
    colnames(distMatrix) <- county_names
  }
  # assign county names
  county_distMatrix <- distMatrix %>% rownames()
  
  
  # loop through all dates
  dates <- data[[time_col]] %>% 
    unique() %>% 
    as.character()
  
  moran_list <- list()
  # for each date, compute Moran's I
  for (date in dates[-1]) {
    cases_date <- data[data[[time_col]] == date, ]
    county_df <- cases_date$CountyName
    ordering <- match(county_distMatrix, county_df)
    
    number_cases_date <- cases_date[ordering, ][[cases_col]]
    
    morans_I_values <- array(NA, dim = 100)
    for (r in seq(1, 100)) {
      set.seed(r)
      # permutate case numbers 
      permutate <- sample(seq(1, 26), 
                          size = 26, 
                          replace = FALSE) 
      
      # compute Moran's I for permutated cases 
      morans_I_values[r] <- ape::Moran.I(number_cases_date[permutate], 
                                         distMatrix, 
                                         scaled = FALSE, 
                                         na.rm = FALSE,
                                         alternative = "two.sided")$observed
    }
    
    quantile_morans_I <- quantile(morans_I_values, 
                                  probs = c(0.025, 0.5, 0.975)) %>% 
      unname()
    
    orig_result <- ape::Moran.I(number_cases_date, 
                                distMatrix, 
                                scaled = FALSE, 
                                na.rm = FALSE,
                                alternative = "two.sided") %>% 
      list.cbind()
    
    # compute Moran's I similar metric based on Kendall's correlation coefficient 
    # for (county in county_distMatrix) {
    #   number_cases_date_county <- number_cases_date[]
    #   
    #   
    # }
    
    moran_list[[date]] <- cbind(orig_result, 
                                "lower_ci" = quantile_morans_I[1], 
                                "upper_ci" = quantile_morans_I[3], 
                                "median" = quantile_morans_I[2])
  }
  
  
  moran_df <- moran_list %>% list.rbind() %>% as.data.frame()
  moran_df$dates <- dates[-1] %>% as.Date()
  
  # p-value 
  morans_p <- moran_df %>% 
    mutate(p = ifelse(observed > upper_ci | observed < lower_ci, 1, 0)) %>% 
    pull(p) %>% 
    mean()
  
  # Visualize and save plot 
  ggplot(moran_df, 
         aes(x = dates, 
             y = observed)) +
    geom_line() +
    xlab("Time") +
    ylab("Moran's I") +
    geom_vline(aes(xintercept = as.Date("27.02.2020",
                                        format = "%d.%m.%Y"), 
                   color = "Initial lockdown")) +
    geom_vline(aes(xintercept = as.Date("18.08.2020",
                                        format = "%d.%m.%Y"), 
                   color = "County-specific restrictions")) +
    geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                        format = "%d.%m.%Y"), 
                   color = "Level-5 lockdown")) +
    geom_vline(aes(xintercept = as.Date("10.05.2021",
                                        format = "%d.%m.%Y"), 
                   color = "Inter-County travel")) +
    geom_vline(aes(xintercept = as.Date("06.03.2022",
                                        format = "%d.%m.%Y"), 
                   color = "End")) +
    geom_line(aes(x = dates, 
                  y = lower_ci), 
              linetype = "dashed", color = "#3A3B3C") +
    geom_line(aes(x = dates, 
                  y = upper_ci), 
              linetype = "dashed", color = "#3A3B3C") +
    geom_line(aes(x = dates, 
                  y = median), 
              linetype = "dashed", color = "#3A3B3C") +
    scale_color_brewer(palette = "Set1") + 
    theme(legend.position = "None")
  ggsave(paste0("Figures/MoransI/covid_moran_", name, ".pdf", collapse = ""), 
         width = 27, height = 14, unit = "cm")
  
  return(morans_p)
}


# Stationarity ------------------------------------------------------------
compute_box_cox <- function(i, 
                            data = datasets_list) {
  data <- data[[i]] %>% 
    as.data.frame() %>% 
    mutate(time = rownames(data[[i]])) %>% 
    gather("CountyName", "COVID_ID", -time)
  
  required_shift <- data %>% pull(COVID_ID) %>% min() %>% abs() + 1
  
  data$COVID_ID <- data$COVID_ID + required_shift
  return(data)
}

# Great Circle distance  --------------------------------------------------
# compute Great Circle distance between all counties 
# requires data as n x 2 matrix with rownames, columns are longitude and 
# latitude coordinates 
circle_distance <- function(data) { 
  pairwise_dist <- matrix(nrow = 26, ncol = 26)
  for (i in seq(1, dim(data)[1])) {
    for (j in seq(1, dim(data)[1])) {
      long1 <- data[i, 1]
      long2 <- data[j, 1]
      lat1 <- data[i, 2]
      lat2 <- data[j, 2]
      pairwise_dist[i, j] <- distm(x = c(long1, lat1), 
                                   y = c(long2, lat2), 
                                   fun = distHaversine)
    }
  }
  # transform into data frame 
  dist_df <- pairwise_dist %>% as.data.frame()
  rownames(dist_df) <- rownames(data)
  colnames(dist_df) <- rownames(data)
  
  # substitute zero diagonals with NA  
  dist_df[dist_df == 0] <- NA
  
  return(dist_df)
}

# ARIMA -------------------------------------------------------------------
# fit ARIMA model, predict the last weeks according to the forecasting 
# window and compute MASE
fit_and_predict_arima <- function(counties = c("Dublin",
                                               "Wicklow", 
                                               "Kerry", 
                                               "Donegal"), 
                                  forecast_window = 10, 
                                  results = results_arima, 
                                  data = covid_cases_df) {
  
  arima_mase <- data.frame("time" = as.Date(NA), 
                           "CountyName" = NA, 
                           "true" = NA, 
                           "predicted" = NA, 
                           "res" = NA, 
                           "mase" = NA, 
                           "type" = NA)
  
  start_date_year <- data %>% 
    rownames() %>% 
    as.Date() %>% 
    min() %>% 
    substr(start = 1, stop = 4) %>% 
    as.numeric()
  
  start_date_month <- data %>% 
    rownames() %>% 
    as.Date() %>% 
    min() %>% 
    substr(start = 6, stop = 7) %>% 
    as.numeric()
  
  l <- data %>% nrow()
  
  for (county in counties) {
    covid_cases_county <- data %>% 
      dplyr::select(county %>% all_of()) %>% 
      ts(frequency = 52, 
         start = c(start_date_year, start_date_month))
    
    # fit ARIMA model, according to best fit model  
    model_county <- results[[county]][[1]]$arma[c(1, 3, 2)]
    
    arima_mod <- arima(covid_cases_county[-((l - forecast_window + 1) : l), ],
                       order = model_county)
    
    # predict based on ARIMA model 
    prediction <- data.frame(time = rownames(data)[((l - forecast_window + 1) : l)], 
                             CountyName = county, 
                             true = data %>% 
                               dplyr::pull(county %>% all_of()) %>% 
                               tail(forecast_window), 
                             predicted = predict(arima_mod, 
                                                 n.ahead = forecast_window)$pred, 
                             mase = 0, 
                             type = "ARIMA") %>% 
      mutate(res = true - predicted)
    prediction$res <- prediction$true - prediction$predicted
    
    
    # compute denominator for MASE 
    denominator <- diff(prediction$true, lag = 1) %>% 
      abs() %>% 
      mean()
    
    # compute MASE 
    for (i in seq(1, nrow(prediction))) {
      prediction[i, ]$mase <- abs(prediction[i, ]$res) / denominator
    }
    arima_mase <- rbind(arima_mase, prediction)
  }
  return(arima_mase %>% na.omit())
}


# GNAR models --------------------------------------------------------------
# fit GNAR model with alpha and beta order according to input arguments  
fit_and_predict <- function(alpha, beta, 
                            globalalpha, 
                            net,
                            vts = covid_cases, 
                            numeric_vertices = FALSE,
                            
                            # if not NULL, coefficients are computed for 
                            # vertex classes  
                            weight_factor = NULL, 
                            
                            # if TRUE, fits INV-D weighting 
                            inverse_distance = TRUE, 
                            
                            # Great circle distance matrix with county names as row- / colnames 
                            distance_matrix = dist_urbanisation %>% as.matrix(), 
                            
                            # data frame with column CountyName and column weight
                            # for population_weight, computes PB weighting 
                            weight_index = population_weight, 
                            
                            # data frame with column CountyName and its numerical encoding 
                            county_index = NULL, 
                            
                            # if TRUE, the original GNARfit() function is applied
                            old = FALSE, 
                            
                            return_model = FALSE, 
                            forecast_window = 10
) {
  
  train_window <- dim(vts)[1] - forecast_window
  
  # fit model according to given settings 
  if (weight_factor %>% is.null()) {
    if (!old) {
      model <- GNARfit_weighting(vts = vts[1:train_window, ], 
                                 net = net, 
                                 numeric_vertices = numeric_vertices, 
                                 alphaOrder = alpha, 
                                 betaOrder = beta, 
                                 globalalpha = globalalpha, 
                                 inverse_distance = inverse_distance,
                                 distance_matrix = distance_matrix,
                                 weight_index = weight_index, 
                                 county_index = county_index
      )
    } 
    if (old) {
      model <- GNARfit(vts = vts[1:train_window, ], 
                       net = net,
                       alphaOrder = alpha, 
                       betaOrder = beta, 
                       globalalpha = globalalpha
      )
    }
    
  } else {
    if (!old) {
      model <- GNARfit_weighting(vts = vts[1:train_window, ],
                                 net = net, 
                                 numeric_vertices = numeric_vertices, 
                                 alphaOrder = alpha, 
                                 betaOrder = beta, 
                                 globalalpha = globalalpha, 
                                 fact.var = weight_factor, 
                                 inverse_distance = inverse_distance,
                                 distance_matrix = distance_matrix,
                                 weight_index = weight_index, 
                                 county_index = county_index
      )
    }
    if (old) {
      model <- GNARfit(vts = vts[1:train_window, ], 
                       net = net,
                       alphaOrder = alpha, 
                       betaOrder = beta, 
                       globalalpha = globalalpha, 
                       fact.var = weight_factor
      )
    }
  }
  
  if (!return_model) {
    # return data frame with RSS and BIC value for model 
    return(data.frame("RSS" = model$mod$residuals^2 %>% sum(), 
                      "BIC" = BIC(model)))
  } 
  if (return_model) {
    # return model 
    return(model)
  }
}

# circle through different alpha and beta order combinations and fit 
# GNAR models with the function "fit_and_predict" 
fit_and_predict_for_many <- function(alpha_options = seq(1, 5), 
                                     beta_options = list(0, 1, 
                                                         c(1, 0), 
                                                         c(1, 1), 
                                                         c(2, 0), 
                                                         c(2, 1), 
                                                         c(2, 2), 
                                                         c(1, 0, 0), 
                                                         c(1, 1, 0),
                                                         c(1, 1, 1), 
                                                         c(2, 1, 1),
                                                         c(1, 0, 0, 0), 
                                                         c(1, 1, 0, 0), 
                                                         c(1, 1, 1, 0), 
                                                         c(1, 1, 1, 1), 
                                                         c(2, 1, 1, 1),
                                                         c(2, 2, 1, 1),
                                                         c(2, 2, 2, 1), 
                                                         c(1, 0, 0, 0, 0), 
                                                         c(1, 1, 0, 0, 0), 
                                                         c(1, 1, 1, 0, 0), 
                                                         c(1, 1, 1, 1, 0), 
                                                         c(1, 1, 1, 1, 1),
                                                         c(2, 1, 1, 1, 1), 
                                                         c(2, 2, 1, 1, 1), 
                                                         c(2, 2, 2, 1, 1)
                                     ), # in form of list
                                     globalalpha = c("TRUE", "FALSE"), 
                                     net, vts = covid_cases,  
                                     numeric_vertices = FALSE, 
                                     
                                     # if not NULL, coefficients are computed for 
                                     # vertex classes
                                     weight_factor = NULL, 
                                     
                                     # if TRUE, fits INV-D weighting 
                                     inverse_distance = TRUE, 
                                     
                                     # Great circle distance matrix with county names as row- / colnames 
                                     distance_matrix = dist_urbanisation %>% as.matrix(),
                                     
                                     # data frame with column CountyName and column weight
                                     # if population_weight, fits PB weighting 
                                     weight_index = population_weight, 
                                     
                                     # data frame with column CountyName and its numerical encoding 
                                     county_index = NULL, 
                                     
                                     # if TRUE, the original GNARfit() function is applied
                                     old = FALSE,
                                     
                                     forecast_window = 10) {
  
  train_window <- dim(vts)[1] - forecast_window
  
  # create all possible combinations of alpha and beta order 
  model_options <-  expand.grid(alpha_options, 
                                beta_options, 
                                globalalpha)
  
  # filter out valid parameter combinations 
  model_options$valid <- (model_options$Var1 == model_options$Var2 %>% 
                            lapply(FUN = function(i) {length(i)}) %>% 
                            unlist())
  model_options_valid <- model_options %>% filter(valid == TRUE)
  
  # for fully connected graph, only the first stage neighbourhood can be 
  # constructed 
  if (any(net %>% 
          GNARtoigraph() %>% 
          igraph::degree() == net$edges %>% 
          length() - 1)) {
    only_1 <- model_options_valid %>% 
      pull(Var2) %>% 
      lapply(FUN = function(i) {not(2 %in% i)}) %>% 
      unlist() 
    
    model_options_valid <- model_options_valid[only_1, ]
  }
  
  BIC_RSS <- data.frame()
  
  for (i in seq(1, nrow(model_options_valid))) {
    # select model settings 
    model_setting <- model_options_valid[i, ]
    
    # create name 
    name <- paste("GNAR", 
                  model_setting$Var1 %>% as.character(), 
                  model_setting$Var2[[1]] %>% paste0(collapse = ""), 
                  model_setting$Var3 %>% as.character(), 
                  sep = "-")
    
    # fit model
    results <- fit_and_predict(alpha = model_setting$Var1, 
                               beta = model_setting$Var2[[1]], 
                               globalalpha = model_setting$Var3 %>% as.logical(), 
                               net = net, vts = vts, 
                               numeric_vertices = numeric_vertices, 
                               weight_factor = weight_factor, 
                               inverse_distance = inverse_distance,
                               distance_matrix = distance_matrix,
                               weight_index = weight_index, 
                               county_index = county_index, 
                               old = old, 
                               forecast_window = forecast_window)
    results$name <- name
    BIC_RSS <- rbind(BIC_RSS, 
                     results)
    
  }
  
  return(BIC_RSS[, c(3, 1, 2)])
}


# fit GNAR models according to certain beta- and alpha-order for all data subsets
fit_subset_models <- function(net, 
                              beta_list, 
                              alpha_vector) {
  
  mod_list <- list()
  
  for (i in seq(1, 5)) {
    mod <- fit_and_predict(alpha = alpha_vector[i], 
                           beta = beta_list[[i]], 
                           net = net, 
                           vts = datasets_list[[i]], 
                           globalalpha = TRUE, 
                           old = TRUE,
                           forecast_window = 0, 
                           return_model = TRUE)
    mod_list[[i]] <- mod
  } 
  return(mod_list)
}


# return model with lowest BIC from a generated results list 
return_best_model <- function(results_list) {
  return(results_list[which.min(results_list$BIC), c(1, 3)])
}

# return list of best performing models for KNN / DNN
return_best_knn_dnn <- function(df) {
  
  res_id <- df[which.min(df$best.model.id.BIC), ][, c(1, 2, 3)]
  res_pop <-df[which.min(df$best.model.pop.BIC), ][, c(1, 4, 5)]
  res_id_class <- df[which.min(df$best.model.id.class.BIC), ][, c(1, 6, 7)]
  res_pop_class <- df[which.min(df$best.model.pop.class.BIC), ][, c(1, 8, 9)]
  res_old <- df[which.min(df$best.model.old.BIC), ][, c(1, 10, 11)]
  res_old_class <- df[which.min(df$best.model.old.class.BIC), ][, c(1, 12, 13)]
  
  return(mapply(c, 
                res_id, 
                res_id_class,
                res_pop,  
                res_pop_class, 
                res_old, 
                res_old_class
  ) %>% data.frame())
}

# Residual analysis -------------------------------------------------------
# compute and plot residuals for GNAR model 
check_and_plot_residuals <- function(model, 
                                     network_name, 
                                     alpha = 4, 
                                     n_ahead = 10, 
                                     counties = c("Dublin",
                                                       "Wicklow", 
                                                       "Kerry", 
                                                       "Donegal"), 
                                     data = covid_cases) {
  
  fitted_df <- model %>% residuals() %>% data.frame()
  colnames(fitted_df) <- colnames(data)
  
  N <- nrow(data)
  # assign time column 
  fitted_df$time <- rownames(data)[-c(1 : alpha, (N - n_ahead + 1) : N)]
  
  fitted_df <- fitted_df %>% gather("CountyName", "residuals", -time)
  
  for (county in counties) {
    # residual QQ plot for selected counties
    g6 <- ggplot(data = fitted_df %>% 
                   filter(CountyName == county),
                 aes(sample = residuals)) +
      geom_qq() +
      geom_qq_line() +
      xlab("theor. quantiles") +
      ylab("emp. quantiles")
    ggsave(filename = paste0("Figures/GNAR_entire_dataset/qq_", 
                             network_name, "_county_",  
                             county, ".pdf"), 
           plot = g6, 
           width = 10, height = 10, unit = "cm")
  }
  return(fitted_df)
}


# compute and plot residuals for GNAR model 
check_and_plot_residuals_subset <- function(model, 
                                     network_name, 
                                     alpha = 5, 
                                     n_ahead = 5, 
                                     counties = c("Dublin",
                                                       "Wicklow", 
                                                       "Kerry", 
                                                       "Donegal"), 
                                     data = covid_cases, 
                                     dataset_name = "GNAR_subsets") {
  
  N <- nrow(data)
  
  residual_list <- list()
  for (i in seq(0, N - n_ahead  - 1)) {
    fitted <- predict(model, newdata = data[1:(n_ahead + i), ])
    colnames(fitted) <- colnames(data)
    real <- data[n_ahead + i + 1, ]
    
    residual_list[[i + 1]] <- real - fitted
  }
  fitted_df <- do.call(rbind.data.frame, residual_list)
  
  # assign time column 
  fitted_df$time <- rownames(data)[-c((N - n_ahead + 1) : N)]
  
  fitted_df <- fitted_df %>% gather("CountyName", "residuals", -time)
  
  for (county in counties) {
  g6 <- ggplot(data = fitted_df %>% 
                   filter(CountyName == county),
                 aes(sample = residuals)) +
      geom_qq() +
      geom_qq_line() +
      xlab("theor. quantiles") +
      ylab("emp. quantiles")
    ggsave(filename = paste0("Figures/",
                             dataset_name, 
                             "/qq_", 
                             network_name, "_county_",  
                             county, 
                             ".pdf"), 
           plot = g6, 
           width = 10, height = 10, unit = "cm")
  }
  return(fitted_df)
}

# KS test for normality of residuals 
ks_residuals <- function(res) {
  res %>% 
    split(res$CountyName) %>% 
    lapply(FUN = function(i) {
      return(ks.test(i$res, "pnorm")$p.value <= 0.025)
    }) %>% 
    list.cbind() %>% 
    table() %>% 
    return()
}


# compute autocorrelation in residuals for each county 
autocorrelation <- function(res, 
                            df_alpha = 5) {
  
  results <- matrix(NA, nrow = 26, ncol = 3)
  i <- 1
  for (county in res$CountyName %>% unique()) {
    
    res_county <- res %>% 
      filter(CountyName == county)
    
    lag <- 1
    
    significant <- TRUE
    while (significant) {
      BL_test <- Box.test(x = res_county$res, 
                          type = "Ljung-Box",
                          lag = df_alpha + lag, 
                          fitdf = df_alpha)
      
      if (df_alpha + lag >= length(res_county$res)) {
        results[i, ] <- c(county, 
                          NA, 
                          NA)
        
        break()
      }
      
      if (BL_test$p.value > 0.05) {
        significant <- FALSE
        
        results[i, ] <- c(county, 
                          df_alpha + lag, 
                          BL_test$p.value)
        break()
      }
      lag <- lag + 1
    }
    i <- i + 1
  }
  return(results)
}


# Restrictions - best model -----------------------------------------------
# analyse the change in coefficient values for best performing GNAR models 
# if fitted to different data subsets 
parameter_development <- function(data_list = datasets_list, 
                                  net, 
                                  network_name,
                                  numeric_vertices = FALSE, 
                                  county_index, 
                                  old = TRUE, 
                                  inverse_distance = FALSE, 
                                  forecast_window = 0, 
                                  alpha = 4, 
                                  beta = c(2, 2, 1, 1), 
                                  globalalpha = TRUE,
                                  return_model = TRUE) {
  param <- list()
  type_of_restriction <- c("Start",
                           "County-sp. lockdowns", 
                           "Level-5 lockdown", 
                           "Inter-county travel", 
                           "End")
  
  # circle through all data subsets and fit the best performing GNAR model
  for (i in seq(1, length(data_list))) {
    mod <- fit_and_predict(vts = data_list[[i]], 
                           net = net, 
                           county_index = county_index, 
                           old = old, 
                           inverse_distance = inverse_distance, 
                           forecast_window = forecast_window, 
                           alpha = alpha, 
                           beta = beta, 
                           globalalpha = globalalpha,
                           return_model = return_model, 
                           numeric_vertices = numeric_vertices)
    
    info <- data.frame("estimate" = mod %>% coef(), 
                       "coefficient" = mod %>% coef() %>% names(), 
                       "restriction" = factor(type_of_restriction[i], 
                                              levels = c("Start",
                                                         "County-sp. lockdowns", 
                                                         "Level-5 lockdown", 
                                                         "Inter-county travel", 
                                                         "End")))
    
    param[[i]] <- cbind(info, 
                        "lower" = confint(mod)[, 1], 
                        "upper" = confint(mod)[, 2]
                        )
    
  }
  
  param_df <- param %>% list.rbind() %>% data.frame()
  
  # visualise the change in coefficient values across data subsets  
  g_alpha <- ggplot(param_df  %>% 
                      filter(grepl("alpha", coefficient)), 
                    aes(x = restriction, 
                        y = estimate, 
                        group = coefficient, 
                        color = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient, 
                      color = coefficient), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype  = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("Figures/ParameterDevelopment/alpha_order_", network_name, ".pdf"), 
         plot = g_alpha, 
         width = 26, height = 14, unit = "cm")
  
  
  g_beta <- ggplot(param_df  %>% 
                     filter(grepl("beta", coefficient)), 
                   aes(x = restriction, 
                       y = estimate, 
                       group = coefficient, 
                       color = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient, 
                      color = coefficient), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("Figures/ParameterDevelopment/beta_order_", network_name, ".pdf"), 
         plot = g_beta, 
         width = 26, height = 14, unit = "cm")
  
  return(list("alpha" = g_alpha, 
              "beta" = g_beta, 
              "dataset" = param_df))
  
}


parameter_development_subsets <- function(data_list = datasets_list, 
                                          net_list, 
                                          numeric_vertices = FALSE, 
                                          county_index, 
                                          old = TRUE, 
                                          inverse_distance = FALSE, 
                                          forecast_window = 0, 
                                          alpha_vector, 
                                          beta_list, 
                                          globalalpha = TRUE,
                                          return_model = TRUE) {
  param <- list()
  type_of_restriction <- c("Start",
                           "County-sp. lockdowns", 
                           "Level-5 lockdown", 
                           "Inter-county travel", 
                           "End")
  
  # circle through all data subsets and fit the best performing GNAR model
  for (i in seq(1, length(data_list))) {
    mod <- fit_and_predict(vts = data_list[[i]], 
                           net = net_list[[i]], 
                           county_index = county_index, 
                           old = old, 
                           inverse_distance = inverse_distance, 
                           forecast_window = forecast_window, 
                           alpha = alpha_vector[i], 
                           beta = beta_list[[i]], 
                           globalalpha = globalalpha,
                           return_model = return_model, 
                           numeric_vertices = numeric_vertices)
    
    info <- data.frame("estimate" = mod %>% coef(), 
                       "coefficient" = mod %>% coef() %>% names(), 
                       "restriction" = factor(type_of_restriction[i], 
                                              levels = c("Start",
                                                         "County-sp. lockdowns", 
                                                         "Level-5 lockdown", 
                                                         "Inter-county travel", 
                                                         "End")))
    param[[i]] <- cbind(info, 
                        "lower" = confint(mod)[, 1], 
                        "upper" = confint(mod)[, 2]
    )
    
  }
  
  param_df <- param %>% list.rbind() %>% data.frame()
  
  # visualise the change in coefficient values across data subsets  
  g_alpha <- ggplot(param_df  %>% 
                      filter(grepl("alpha", coefficient)), 
                    aes(x = restriction, 
                        y = estimate, 
                        group = coefficient, 
                        color = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient, 
                      color = coefficient), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype  = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = "Figures/ParameterDevelopment/alpha_order_subsets.pdf", 
         plot = g_alpha, 
         width = 26, height = 14, unit = "cm")
  
  
  g_beta <- ggplot(param_df  %>% 
                     filter(grepl("beta", coefficient)), 
                   aes(x = restriction, 
                       y = estimate, 
                       group = coefficient, 
                       color = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient, 
                      color = coefficient), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = "Figures/ParameterDevelopment/beta_order_subsets.pdf", 
         plot = g_beta, 
         width = 26, height = 14, unit = "cm")
  
  return(list("alpha" = g_alpha, 
              "beta" = g_beta, 
              "dataset" = param_df))
  
}

parameter_development_phases <- function(data_list = datasets_list_coarse, 
                                          net_list, 
                                          numeric_vertices = FALSE, 
                                          county_index, 
                                          old = TRUE, 
                                          inverse_distance = FALSE, 
                                          forecast_window = 0, 
                                          alpha_vector, 
                                          beta_list, 
                                          globalalpha = TRUE,
                                          return_model = TRUE, 
                                         name = "phases") {
  param <- list()
  type_of_restriction <- c("Restricted",
                           "Unrestricted")
  
  # circle through all data subsets and fit the best performing GNAR model
  for (i in seq(1, length(data_list))) {
    mod <- fit_and_predict(vts = data_list[[i]], 
                           net = net_list[[i]], 
                           county_index = county_index, 
                           old = old, 
                           inverse_distance = inverse_distance, 
                           forecast_window = forecast_window, 
                           alpha = alpha_vector[i], 
                           beta = beta_list[[i]], 
                           globalalpha = globalalpha,
                           return_model = return_model, 
                           numeric_vertices = numeric_vertices)
    
    info <- data.frame("estimate" = mod %>% coef(), 
                       "coefficient" = mod %>% coef() %>% names(), 
                       "restriction" = factor(type_of_restriction[i], 
                                              levels = c("Restricted",
                                                         "Unrestricted")))
    param[[i]] <- cbind(info, 
                        "lower" = confint(mod)[, 1], 
                        "upper" = confint(mod)[, 2]
    )
    
  }
  
  param_df <- param %>% 
    list.rbind() %>% 
    data.frame()
  
  # format coefficient names 
  param_df$coefficient_formated <- gsub("dmat2", "", 
                                        param_df$coefficient)
  
  # visualise the change in coefficient values across data subsets  
  g_alpha <- ggplot(param_df  %>% 
                      filter(grepl("alpha", coefficient_formated)), 
                    aes(x = restriction, 
                        y = estimate, 
                        group = coefficient_formated, 
                        color = coefficient_formated)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient_formated, 
                      color = coefficient_formated), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype  = "dashed") +
    labs(x = "Pandemic phase", 
         y = "coefficient values", 
         color = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") 
  # +
  #   guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("Figures/ParameterDevelopment/alpha_order_", 
                       name,
                       ".pdf", 
                       collapse = ""), 
         plot = g_alpha, 
         width = 18, height = 15, unit = "cm")
  
  
  g_beta <- ggplot(param_df  %>% 
                     filter(grepl("beta", coefficient_formated)), 
                   aes(x = restriction, 
                       y = estimate, 
                       group = coefficient_formated, 
                       color = coefficient_formated)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient_formated, 
                      color = coefficient_formated), 
                  width=0.2, 
                  alpha = 0.5) +
    geom_line(linetype = "dashed") +
    labs(x = "Pandemic phase", 
         y = "coefficient values", 
         color = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") 
  # +
  #   guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("Figures/ParameterDevelopment/beta_order_", 
                       name,
                       ".pdf", 
                       collapse = ""), 
         plot = g_beta, 
         width = 18, height = 15, unit = "cm")
  
  return(list("alpha" = g_alpha, 
              "beta" = g_beta, 
              "dataset" = param_df))
  
}


# Restrictions - tuning parameters ----------------------------------------
# circle through different alpha and beta order combinations and fit 
# GNAR models with the function "fit_and_predict" for each data subset 
fit_and_predict_for_restrictions <- function(alpha_options = seq(1, 5), 
                                             beta_options = list(0, 1, 
                                                                 c(1, 0), 
                                                                 c(1, 1), 
                                                                 c(2, 0), 
                                                                 c(2, 1), 
                                                                 c(2, 2), 
                                                                 c(1, 0, 0), 
                                                                 c(2, 0, 0),
                                                                 c(3, 0, 0),
                                                                 c(4, 0, 0),
                                                                 c(5, 0, 0),
                                                                 c(1, 1, 0),
                                                                 c(2, 1, 0),
                                                                 c(3, 1, 0),
                                                                 c(4, 1, 0),
                                                                 c(5, 1, 0),
                                                                 c(2, 2, 0),
                                                                 c(1, 1, 1), 
                                                                 c(2, 1, 1),
                                                                 c(3, 1, 1),
                                                                 c(4, 1, 1),
                                                                 c(5, 1, 1),
                                                                 c(2, 2, 1), 
                                                                 c(2, 2, 2),
                                                                 c(3, 2, 2), 
                                                                 c(4, 2, 2),
                                                                 c(5, 2, 2), 
                                                                 c(1, 0, 0, 0), 
                                                                 c(2, 0, 0, 0),
                                                                 c(3, 0, 0, 0),
                                                                 c(4, 0, 0, 0),
                                                                 c(5, 0, 0, 0),
                                                                 c(1, 1, 0, 0), 
                                                                 c(2, 1, 0, 0),
                                                                 c(3, 1, 0, 0),
                                                                 c(4, 1, 0, 0),
                                                                 c(5, 1, 0, 0),
                                                                 c(2, 2, 0, 0),
                                                                 c(1, 1, 1, 0), 
                                                                 c(2, 1, 1, 0),
                                                                 c(3, 1, 1, 0), 
                                                                 c(4, 1, 1, 0),
                                                                 c(5, 1, 1, 0),
                                                                 c(2, 2, 1, 0), 
                                                                 c(2, 2, 2, 0),
                                                                 c(3, 2, 2, 0),
                                                                 c(4, 2, 2, 0),
                                                                 c(5, 2, 2, 0),
                                                                 c(1, 1, 1, 1), 
                                                                 c(2, 1, 1, 1),
                                                                 c(3, 1, 1, 1), 
                                                                 c(4, 1, 1, 1),
                                                                 c(5, 1, 1, 1),
                                                                 c(2, 2, 1, 1),
                                                                 c(2, 2, 2, 1), 
                                                                 c(3, 2, 2, 1), 
                                                                 c(4, 2, 2, 1),
                                                                 c(5, 2, 2, 1),
                                                                 c(2, 2, 2, 2),
                                                                 c(1, 0, 0, 0, 0),
                                                                 c(2, 0, 0, 0, 0),
                                                                 c(3, 0, 0, 0, 0), 
                                                                 c(4, 0, 0, 0, 0),
                                                                 c(5, 0, 0, 0, 0), 
                                                                 c(6, 0, 0, 0, 0), 
                                                                 c(7, 0, 0, 0, 0), 
                                                                 c(1, 1, 0, 0, 0),
                                                                 c(2, 1, 0, 0, 0), 
                                                                 c(3, 1, 0, 0, 0),
                                                                 c(4, 1, 0, 0, 0),
                                                                 c(5, 1, 0, 0, 0),
                                                                 c(6, 1, 0, 0, 0),
                                                                 c(7, 1, 0, 0, 0),
                                                                 c(2, 2, 0, 0, 0),
                                                                 c(1, 1, 1, 0, 0), 
                                                                 c(2, 1, 1, 0, 0),
                                                                 c(3, 1, 1, 0, 0), 
                                                                 c(4, 1, 1, 0, 0),
                                                                 c(5, 1, 1, 0, 0),
                                                                 c(6, 1, 1, 0, 0),
                                                                 c(7, 1, 1, 0, 0),
                                                                 c(2, 2, 1, 0, 0), 
                                                                 c(2, 2, 2, 0, 0),
                                                                 c(1, 1, 1, 1, 0), 
                                                                 c(2, 1, 1, 1, 0),
                                                                 c(3, 1, 1, 1, 0), 
                                                                 c(4, 1, 1, 1, 0),
                                                                 c(5, 1, 1, 1, 0),
                                                                 c(6, 1, 1, 1, 0),
                                                                 c(7, 1, 1, 1, 0),
                                                                 c(2, 2, 1, 1, 0), 
                                                                 c(2, 2, 2, 1, 0),
                                                                 c(2, 2, 2, 2, 0),
                                                                 c(3, 2, 2, 2, 0), 
                                                                 c(4, 2, 2, 2, 0),
                                                                 c(5, 2, 2, 2, 0),
                                                                 c(6, 2, 2, 2, 0),
                                                                 c(7, 2, 2, 2, 0),
                                                                 c(1, 1, 1, 1, 1),
                                                                 c(2, 1, 1, 1, 1),
                                                                 c(3, 1, 1, 1, 1), 
                                                                 c(4, 1, 1, 1, 1),
                                                                 c(5, 1, 1, 1, 1),
                                                                 c(6, 1, 1, 1, 1),
                                                                 c(7, 1, 1, 1, 1),
                                                                 c(2, 2, 1, 1, 1), 
                                                                 c(2, 2, 2, 1, 1),
                                                                 c(3, 2, 2, 1, 1),
                                                                 c(4, 2, 2, 1, 1),
                                                                 c(5, 2, 2, 1, 1),
                                                                 c(6, 2, 2, 1, 1),
                                                                 c(7, 2, 2, 1, 1),
                                                                 c(2, 2, 2, 2, 1), 
                                                                 c(2, 2, 2, 2, 2)
                                             ), # in form of list
                                             globalalpha = TRUE, 
                                             net, 
                                             data_list = datasets_list,  
                                             numeric_vertices = FALSE,
                                             
                                             # if not NULL, fit coefficients for 
                                             # vertex classes 
                                             weight_factor = NULL, 
                                             
                                             # if TRUE, INV-D weighting 
                                             inverse_distance = FALSE, 
                                             
                                             # Great circle distance matrix with county names as row- / colnames 
                                             distance_matrix = dist_urbanisation %>% as.matrix(), 
                                             
                                             # data frame with column CountyName and column weight
                                             weight_index = population_weight, 
                                             
                                             # data frame with column CountyName and its numerical encoding 
                                             county_index = NULL, 
                                             
                                             # if TRUE, the original GNARfit() function is applied
                                             old = TRUE,
                                             
                                             forecast_window = 5,
                                             upper_limit = 5) {
  
  # construct valid beta options list 
  exclude <- lapply(beta_options, 
                    FUN = function(i) any(i > upper_limit)) %>% 
    unlist()
  
  beta_valid_options <- beta_options[!exclude]
  
  param <- list()
  
  for (i in seq(1, length(data_list))) {
    res <- fit_and_predict_for_many(vts = data_list[[i]], 
                                    net = net, 
                                    alpha_options = alpha_options, 
                                    beta_options = beta_valid_options, 
                                    globalalpha = globalalpha,
                                    old = old,
                                    forecast_window = forecast_window, 
                                    numeric_vertices = numeric_vertices, 
                                    weight_factor = weight_factor, 
                                    inverse_distance = inverse_distance, 
                                    distance_matrix = distance_matrix, 
                                    weight_index = weight_index,
                                    county_index = county_index)
    
    param[[i]] <- return_best_model(results_list = res)
  }
  param_df <- do.call(rbind.data.frame, param)
  
  param_df$data_subset <- seq(1, length(data_list))
  
  return(param_df[, c(3, 1, 2)])
}


fit_and_predict_for_restrictions_all_models <- function(alpha_options = seq(1, 5), 
                                                        beta_options = list(0, 1, 
                                                                            c(1, 0), 
                                                                            c(1, 1), 
                                                                            c(2, 0), 
                                                                            c(2, 1), 
                                                                            c(2, 2), 
                                                                            c(1, 0, 0), 
                                                                            c(2, 0, 0),
                                                                            c(3, 0, 0),
                                                                            c(4, 0, 0),
                                                                            c(5, 0, 0),
                                                                            c(1, 1, 0),
                                                                            c(2, 1, 0),
                                                                            c(3, 1, 0),
                                                                            c(4, 1, 0),
                                                                            c(5, 1, 0),
                                                                            c(2, 2, 0),
                                                                            c(1, 1, 1), 
                                                                            c(2, 1, 1),
                                                                            c(3, 1, 1),
                                                                            c(4, 1, 1),
                                                                            c(5, 1, 1),
                                                                            c(2, 2, 1), 
                                                                            c(2, 2, 2),
                                                                            c(3, 2, 2), 
                                                                            c(4, 2, 2),
                                                                            c(5, 2, 2), 
                                                                            c(1, 0, 0, 0), 
                                                                            c(2, 0, 0, 0),
                                                                            c(3, 0, 0, 0),
                                                                            c(4, 0, 0, 0),
                                                                            c(5, 0, 0, 0),
                                                                            c(1, 1, 0, 0), 
                                                                            c(2, 1, 0, 0),
                                                                            c(3, 1, 0, 0),
                                                                            c(4, 1, 0, 0),
                                                                            c(5, 1, 0, 0),
                                                                            c(2, 2, 0, 0),
                                                                            c(1, 1, 1, 0), 
                                                                            c(2, 1, 1, 0),
                                                                            c(3, 1, 1, 0), 
                                                                            c(4, 1, 1, 0),
                                                                            c(5, 1, 1, 0),
                                                                            c(2, 2, 1, 0), 
                                                                            c(2, 2, 2, 0),
                                                                            c(3, 2, 2, 0),
                                                                            c(4, 2, 2, 0),
                                                                            c(5, 2, 2, 0),
                                                                            c(1, 1, 1, 1), 
                                                                            c(2, 1, 1, 1),
                                                                            c(3, 1, 1, 1), 
                                                                            c(4, 1, 1, 1),
                                                                            c(5, 1, 1, 1),
                                                                            c(2, 2, 1, 1),
                                                                            c(2, 2, 2, 1), 
                                                                            c(3, 2, 2, 1), 
                                                                            c(4, 2, 2, 1),
                                                                            c(5, 2, 2, 1),
                                                                            c(2, 2, 2, 2),
                                                                            c(1, 0, 0, 0, 0),
                                                                            c(2, 0, 0, 0, 0),
                                                                            c(3, 0, 0, 0, 0), 
                                                                            c(4, 0, 0, 0, 0),
                                                                            c(5, 0, 0, 0, 0), 
                                                                            c(6, 0, 0, 0, 0), 
                                                                            c(7, 0, 0, 0, 0), 
                                                                            c(1, 1, 0, 0, 0),
                                                                            c(2, 1, 0, 0, 0), 
                                                                            c(3, 1, 0, 0, 0),
                                                                            c(4, 1, 0, 0, 0),
                                                                            c(5, 1, 0, 0, 0),
                                                                            c(6, 1, 0, 0, 0),
                                                                            c(7, 1, 0, 0, 0),
                                                                            c(2, 2, 0, 0, 0),
                                                                            c(1, 1, 1, 0, 0), 
                                                                            c(2, 1, 1, 0, 0),
                                                                            c(3, 1, 1, 0, 0), 
                                                                            c(4, 1, 1, 0, 0),
                                                                            c(5, 1, 1, 0, 0),
                                                                            c(6, 1, 1, 0, 0),
                                                                            c(7, 1, 1, 0, 0),
                                                                            c(2, 2, 1, 0, 0), 
                                                                            c(2, 2, 2, 0, 0),
                                                                            c(1, 1, 1, 1, 0), 
                                                                            c(2, 1, 1, 1, 0),
                                                                            c(3, 1, 1, 1, 0), 
                                                                            c(4, 1, 1, 1, 0),
                                                                            c(5, 1, 1, 1, 0),
                                                                            c(6, 1, 1, 1, 0),
                                                                            c(7, 1, 1, 1, 0),
                                                                            c(2, 2, 1, 1, 0), 
                                                                            c(2, 2, 2, 1, 0),
                                                                            c(2, 2, 2, 2, 0),
                                                                            c(3, 2, 2, 2, 0), 
                                                                            c(4, 2, 2, 2, 0),
                                                                            c(5, 2, 2, 2, 0),
                                                                            c(6, 2, 2, 2, 0),
                                                                            c(7, 2, 2, 2, 0),
                                                                            c(1, 1, 1, 1, 1),
                                                                            c(2, 1, 1, 1, 1),
                                                                            c(3, 1, 1, 1, 1), 
                                                                            c(4, 1, 1, 1, 1),
                                                                            c(5, 1, 1, 1, 1),
                                                                            c(6, 1, 1, 1, 1),
                                                                            c(7, 1, 1, 1, 1),
                                                                            c(2, 2, 1, 1, 1), 
                                                                            c(2, 2, 2, 1, 1),
                                                                            c(3, 2, 2, 1, 1),
                                                                            c(4, 2, 2, 1, 1),
                                                                            c(5, 2, 2, 1, 1),
                                                                            c(6, 2, 2, 1, 1),
                                                                            c(7, 2, 2, 1, 1),
                                                                            c(2, 2, 2, 2, 1), 
                                                                            c(2, 2, 2, 2, 2)
                                                        ), # in form of list
                                                        globalalpha = TRUE, 
                                                        net, 
                                                        data_list = datasets_list,  
                                                        numeric_vertices = FALSE,
                                                        
                                                        # if not NULL, fit coefficients for 
                                                        # vertex classes 
                                                        weight_factor = NULL, 
                                                        
                                                        # if TRUE, INV-D weighting 
                                                        inverse_distance = FALSE, 
                                                        
                                                        # Great circle distance matrix with county names as row- / colnames 
                                                        distance_matrix = dist_urbanisation %>% as.matrix(), 
                                                        
                                                        # data frame with column CountyName and column weight
                                                        weight_index = population_weight, 
                                                        
                                                        # data frame with column CountyName and its numerical encoding 
                                                        county_index = NULL, 
                                                        
                                                        # if TRUE, the original GNARfit() function is applied
                                                        old = TRUE,
                                                        
                                                        forecast_window = 5,
                                                        upper_limit = 5) {
  
  # construct valid beta options list 
  exclude <- lapply(beta_options, 
                    FUN = function(i) any(i > upper_limit)) %>% 
    unlist()
  
  beta_valid_options <- beta_options[!exclude]
  
  param <- list()
  
  for (i in seq(1, length(data_list))) {
    res <- fit_and_predict_for_many(vts = data_list[[i]], 
                                    net = net, 
                                    alpha_options = alpha_options, 
                                    beta_options = beta_valid_options, 
                                    globalalpha = globalalpha,
                                    old = old,
                                    forecast_window = forecast_window, 
                                    numeric_vertices = numeric_vertices, 
                                    weight_factor = weight_factor, 
                                    inverse_distance = inverse_distance, 
                                    distance_matrix = distance_matrix, 
                                    weight_index = weight_index,
                                    county_index = county_index)
    
    param[[i]] <- res %>% 
      mutate(data_subset = rep(i, nrow(res)))
  }
  param_df <- do.call(rbind.data.frame, param)
  
  return(param_df[, c(3, 1, 2, 4)])
}

# MASE --------------------------------------------------------------------
# compute MASE for certain counties 
compute_MASE <- function(model, 
                         network_name, 
                         n_ahead = 10, 
                         counties = c("Dublin", 
                                      "Wicklow", 
                                      "Kerry", 
                                      "Donegal"), 
                         data_df = covid_cases_df) {
  
  predicted_df <- predict(model,  
                          n.ahead = n_ahead) %>% 
    as.data.frame()
  
  colnames(predicted_df) <- colnames(data_df %>% dplyr::select(-time))
  
  length_data_df <- nrow(data_df)
  
  # start date for the period of prediction 
  prediction_time <- data_df[1:(length_data_df - n_ahead), ] %>% 
    rownames() %>% 
    tail(1)
  
  end_date <- data_df[length_data_df, ] %>% 
    rownames()
  
  predicted_df$time <- seq(as.Date(prediction_time) + 7, 
                           as.Date(end_date), 
                           by = 7)
  
  # compare true 1-lag COVID-19 ID and predicted values 
  true <- data_df[(length_data_df - n_ahead + 1):length_data_df, ] %>% 
    gather(key = "CountyName",
           value = "true", 
           -time)
  pred <- predicted_df %>% 
    gather(key = "CountyName",
           value = "predicted", 
           -time)
  # create data frame to assess predictive performance 
  check_predictions_df <- left_join(true, pred, 
                                    by = c("CountyName", "time")) %>% 
    mutate(res = true - predicted)
  
  mase_df <- data.frame("time" = as.Date(NA), 
                        "CountyName" = NA, 
                        "true" = NA, 
                        "predicted" = NA, 
                        "res" = NA, 
                        "mase" = NA, 
                        "type" = NA)
  
  for (county in counties) {
    check_predictions_county <- check_predictions_df %>% 
      filter(CountyName == county) 
    
    check_predictions_county$mase <- 0
    
    # compute denominator for MASE
    denominator <- diff(check_predictions_county$true, lag = 1) %>% 
      abs() %>% 
      mean()
    
    for (i in seq(1, nrow(check_predictions_county))) {
      # compute MASE values 
      check_predictions_county[i, ]$mase <- abs(check_predictions_county[i, ]$res) / denominator
      }
    
    check_predictions_county$type <- network_name
    
    mase_df <- rbind.data.frame(mase_df, 
                                check_predictions_county)
    
  }
  
  return(mase_df %>% na.omit())
}


plot_mase_I <- function(mase_overview, 
                        mase_name = "restricted", 
                        types = c("ARIMA", 
                                  "subset_1_gabriel", 
                                  "subset_1_relative", 
                                  "subset_1_soi", 
                                  "subset_1_delaunay", 
                                  "subset_1_train"),
                        type_name = "delaunay", 
                        counties_subset = all_counties[1:13], 
                        number_counties = 1, 
                        color_types = c("ARIMA" = "#3A3B3C", 
                                        "subset_1_gabriel" = "#00BFC4", 
                                        "subset_1_relative" = "#00B0F6", 
                                        "subset_1_soi" = "#9590FF", 
                                        "subset_1_delaunay" = "#E76BF3", 
                                        "subset_1_train" = "#FF62BC"), 
                        label_networks = c("ARIMA", 
                                           "Gabriel", 
                                           "Relative", 
                                           "SOI", 
                                           "Delaunay", 
                                           "Train"), 
                        lag = "lag_1") {
  # filter current considered network 
  mase_overview_df <- mase_overview %>% 
    filter(type %in% types, 
           CountyName %in% counties_subset) 
  # create factor type 
  mase_overview_df$type <- factor(mase_overview_df$type, 
                                  levels = types)
  
  
  g <- ggplot(mase_overview_df, 
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
    scale_color_manual(values = color_types,
                       labels = label_networks, 
                       name = "Network")
  ggsave(filename = paste0("Figures/GNAR_pandemic_phases/mase_", 
                           type_name,
                           "_", 
                           mase_name, 
                           "_", 
                           lag, 
                           number_counties,
                           ".pdf", 
                           collapse = ""), 
         plot = g,
         width = 30, height = 13, units = "cm")
  
  return(g)
}

plot_mase_II <- function(mase_overview, 
                         mase_name = "restricted", 
                         counties_subset,
                         number_counties, 
                         lag = "lag_1", 
                         type_name = "knn", 
                         types = c("ARIMA", 
                                   "subset_1_knn", 
                                   "subset_1_dnn", 
                                   "subset_1_complete", 
                                   "subset_1_queen", 
                                   "subset_1_eco_hub"),
                         color_types = c("ARIMA" = "#3A3B3C", 
                                         "subset_1_knn" = "#F8766D", 
                                         "subset_1_dnn" = "#D89000", 
                                         "subset_1_complete" = "#A3A500", 
                                         "subset_1_queen" = "#39B600", 
                                         "subset_1_eco_hub" = "#00BF7D")) {
  plot_mase_I(mase_overview = mase_overview, 
              mase_name = mase_name, 
              types = types,
              type_name = type_name, 
              counties_subset = counties_subset, 
              number_counties = number_counties, 
              color_types = color_types, 
              label_networks = c("ARIMA", 
                                 "KNN", 
                                 "DNN", 
                                 "Complete",
                                 "Queen", 
                                 "Eco hub"), 
              lag = lag)
}



# Prediction vs. Fitted ---------------------------------------------------
plot_predicted_vs_fitted_I <- function(mase_overview, 
                                       mase_name = "restricted", 
                                       types = c("ARIMA", 
                                                 "subset_1_gabriel", 
                                                 "subset_1_relative", 
                                                 "subset_1_soi", 
                                                 "subset_1_delaunay", 
                                                 "subset_1_train"),
                                       type_name = "delaunay", 
                                       counties_subset = all_counties[1:13], 
                                       number_counties = 1, 
                                       color_types = c("ARIMA" = "#3A3B3C", 
                                                       "subset_1_gabriel" = "#00BFC4", 
                                                       "subset_1_relative" = "#00B0F6", 
                                                       "subset_1_soi" = "#9590FF", 
                                                       "subset_1_delaunay" = "#E76BF3", 
                                                       "subset_1_train" = "#FF62BC"), 
                                       label_networks = c("ARIMA", 
                                                          "Gabriel", 
                                                          "Relative", 
                                                          "SOI", 
                                                          "Delaunay", 
                                                          "Train"), 
                                       lag = "lag_1") {
  # filter current considered network 
  mase_overview_df <- mase_overview %>% 
    filter(type %in% types, 
           CountyName %in% counties_subset) 
  # create factor type 
  mase_overview_df$type <- factor(mase_overview_df$type, 
                                  levels = types)
  
  
  g <- ggplot(mase_overview_df, 
              aes(x = time, 
                  y = predicted, 
                  color = type)) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_point(aes(x = time, 
                   y = true), 
               color = "darkgrey", 
               shape = 2) +
    geom_line(aes(x = time, 
                  y = true), 
              linetype = "dotted", 
              color = "darkgrey") +
    xlab("Time") + 
    ylab("COVID-19 ID") +
    facet_grid(~ CountyName) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1)) +
    scale_color_manual(values = color_types,
                       labels = label_networks, 
                       name = "Network")
  ggsave(filename = paste0("Figures/GNAR_pandemic_phases/predicted_", 
                           type_name,
                           "_", 
                           mase_name, 
                           "_", 
                           lag, 
                           number_counties,
                           ".pdf", 
                           collapse = ""), 
         plot = g,
         width = 26, height = 13, units = "cm")
  
  return(g)
}

plot_predicted_vs_fitted_II <- function(mase_overview, 
                                        mase_name = "restricted",
                                        counties_subset,
                                        number_counties, 
                                        lag = "lag_1", 
                                        types = c("ARIMA", 
                                                  "subset_1_knn", 
                                                  "subset_1_dnn", 
                                                  "subset_1_complete", 
                                                  "subset_1_queen", 
                                                  "subset_1_eco_hub"),
                                        color_types = c("ARIMA" = "#3A3B3C", 
                                                        "subset_1_knn" = "#F8766D", 
                                                        "subset_1_dnn" = "#D89000", 
                                                        "subset_1_complete" = "#A3A500", 
                                                        "subset_1_queen" = "#39B600", 
                                                        "subset_1_eco_hub" = "#00BF7D")) {
  plot_predicted_vs_fitted_I(mase_overview = mase_overview, 
                             mase_name = mase_name, 
                             types = types,
                             type_name = "knn", 
                             counties_subset = counties_subset, 
                             number_counties = number_counties, 
                             color_types = color_types, 
                             label_networks = c("ARIMA", 
                                                "KNN", 
                                                "DNN", 
                                                "Complete",
                                                "Queen", 
                                                "Eco hub"), 
                             lag = lag)
}


# Simulation data  --------------------------------------------------------
# generate simulated data
simulate_time_series <- function(gnar_object, 
                                 model_coef, # model coefficients
                                 alpha_order,
                                 beta_order, # beta order
                                 initial_mean = 10, 
                                 var_noise, # variance for error term
                                 county_index, 
                                 counties = counties_v, 
                                 timeframe = 120) {
  # determine max stage across lags 
  max_stage <- beta_order %>% max()
  
  # generate white noise 1-lag COVID-19 ID values for first n days
  initial_array <- rnorm(n = alpha_order  * 26, 
                      mean = initial_mean, 
                      sd = sqrt(var_noise)) %>%
    array(dim = c(alpha_order, 26))
  
  # initialize data frame to store simulated data
  simulated_array <- array(0, dim = c(timeframe, 26))
  simulated_array[1:alpha_order, ] <- initial_array
  
  
  # assign counties to columns  
  colnames(simulated_array) <- counties
  
  # current time 
  current_end <- initial_array %>% nrow()
  
  # predict for weeks according to time frame
  while (current_end < timeframe) {
    simulated_values <- array(0, dim = 26)
    
    for (county in counties) {
      df_index <- which(colnames(simulated_array) == county)
      
      # create data frame consisting of the relevant county as column
      initial_df_county <- simulated_array[, county]
      
      # numeric index for vertex according to county index
      county_to_vertex <- county_index[county_index$CountyName == county, ]$index
      
      
      # separate alpha and beta coefficients 
      alpha_coef <- model_coef  %>% 
        filter(grepl("alpha", type))
      
      beta_coef <- model_coef  %>% 
        filter(grepl("beta", type))
      
      # compute shortest path length between county and all remaining counties
      spl_county <- distances(graph = gnar_object %>% 
                                GNARtoigraph(), 
                              v = county_to_vertex) %>% 
        t() %>% 
        as.data.frame()
      colnames(spl_county) <- "spl"
      spl_county$index <- seq(1, 26)
      
      # add shortest path length to county index data frame 
      spl_county_names <- left_join(county_index, 
                                    spl_county, 
                                    by = "index")
      
      # create vector for beta terms to be saved in 
      beta_part <- array(0, dim = c(alpha_order, max_stage))
      
      for (lagged in seq(1, alpha_order)) {
        # determine neighborhood stage for certain lag
        beta_order_component <- beta_order[lagged]
        
        # determine beta coefficients for certain lag 
        beta_coef_lagged <- beta_coef %>%
          filter(grepl(paste0("beta", lagged), type))
        
        # for every stage, identify neighbors and compute term based on beta coefficients and lagged values 
        if (beta_order_component != 0) {
          for (stage in seq(1, beta_order_component)) {
            # identify vertices in certain neighbourhood stage 
            stage_neighbourhood <- spl_county_names %>% 
              filter(spl == stage) %>% 
              pull(CountyName)
            
            # compute SPL weights as the inverse of numbers of vertices
            # in neighbourhood
            neighbour_weight <- 1 / (stage_neighbourhood %>% 
                                       length())
            
            # filter data for neighbours 
            initial_df_neighbours <- simulated_array[, stage_neighbourhood]
            
            # compute beta term for certain time lag 
            beta_part[lagged, stage] <- beta_coef_lagged$param[stage] * neighbour_weight * sum(initial_df_neighbours[current_end + 1 - lagged])
            
          }
        } 
      }
      
      # compute alpha term as sum of all time lags 
      alpha_term <- sum(alpha_coef$param * initial_df_county[(current_end + 1 - alpha_order) : current_end])
      
      # compute over beta term as sum of all time lags 
      beta_term <- beta_part %>% sum()
      
      # compute simulated value plus random error
      simulated_values[df_index] <- alpha_term + beta_term + rnorm(n = 1, 
                                                        mean = 0, 
                                                        sd = sqrt(var_noise))
      
    }
    # add new row with simulated data points for each county  
    simulated_array[current_end + 1, ] <- simulated_values
    
    
    # compute current time end point 
    current_end <- current_end + 1
  }
  
  # assign fictional "time" 
  simulated_df <- simulated_array %>% 
    as.data.frame() %>% 
    mutate(time = seq(1, timeframe))
  
  # transform into long data table 
  simulated_df_long <- simulated_df %>% 
    gather("CountyName", 
           "ID", 
           -time)
  
  return(simulated_df_long)
}


# re-compute GNAR model coefficients for simulated data 
reconstruct_coefficients <- function(model_coef,
                                     simulation_df, 
                                     alphaOrder, 
                                     betaOrder, 
                                     gnar_object) {

  simulation_short <- simulation_df %>% 
    spread(CountyName, ID) %>% 
    dplyr::select(-time) %>% 
    as.matrix()
  
  gnar_model <- GNARfit(vts = simulation_short, 
                        net = gnar_object, 
                        alphaOrder = alphaOrder, 
                        betaOrder = betaOrder, 
                        globalalpha = TRUE)
  
  # extract estimated model coefficients 
  fitted_coef_df <- gnar_model %>% 
    coef() %>% 
    as.data.frame()
  colnames(fitted_coef_df) <- "recomp_param"
  
  fitted_coef_df <- fitted_coef_df %>% 
    mutate(type = rownames(fitted_coef_df))
  
  # join true and estimated coefficients 
  compare_coef <- left_join(fitted_coef_df, 
                            model_coef, 
                            by = "type")
    
  
  # compute 95% confidence interval 
  ci_coef <- gnar_model$mod %>% confint() 
  
  # add confidence interval to coefficient estimate 
  compare_coef$recomp_param_ci <- paste0(compare_coef$recomp_param %>% round(2), 
                                    " [", 
                                    ci_coef[, 1] %>% round(2),
                                    ", ",
                                    ci_coef[, 2] %>% round(2),
                                    "]")
  return(compare_coef[, c(2, 3, 4)])
}



