### FUNCTIONS

# source adapted functions from GNAR package 
source("GNAR/GNARdesign.R")
source("GNAR/GNARfit.R")
source("GNAR/NofNeighbours.R")

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

# Knox Mantel test 
rank_moran_I_permutation_test <- function(data = COVID_weekly_data,
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
  
  spearman_corr_list <- list()
  # for each date, compute Moran's I
  for (date in dates[-1]) {
    cases_date <- data[data[[time_col]] == date, ]
    county_df <- cases_date$CountyName
    ordering <- match(county_distMatrix, county_df)
    
    
    rank_cases_date <- cases_date[ordering, ][[cases_col]] %>% 
      rank()
    names(rank_cases_date) <- county_distMatrix
    
    spearman_corr_values <- array(NA, dim = 100)
    for (r in seq(1, 100)) {
      set.seed(r)
      # permutate case numbers 
      permutate <- sample(seq(1, 26), 
                          size = 26, 
                          replace = FALSE) 
      
      # compute Moran's I for permutated cases 
      spearman_corr_values[r] <- ape::Moran.I(rank_cases_date[permutate], 
                                         distMatrix, 
                                         scaled = FALSE, 
                                         na.rm = FALSE,
                                         alternative = "two.sided")$observed
    }
    
    quantile_spearman_corr <- quantile(spearman_corr_values, 
                                  probs = c(0.025, 0.5, 0.975)) %>% 
      unname()
    
    orig_result <- ape::Moran.I(rank_cases_date, 
                                distMatrix, 
                                scaled = FALSE, 
                                na.rm = FALSE,
                                alternative = "two.sided") %>% 
      list.cbind()
    
    spearman_corr_list[[date]] <- cbind(orig_result, 
                                "lower_ci" = quantile_spearman_corr[1], 
                                "upper_ci" = quantile_spearman_corr[3], 
                                "median" = quantile_spearman_corr[2])
  }
  
  
  spearman_corr_df <- spearman_corr_list %>% list.rbind() %>% as.data.frame()
  spearman_corr_df$dates <- dates[-1] %>% as.Date()
  
  # p-value 
  spearman_corr_p <- spearman_corr_df %>% 
    mutate(p = ifelse(observed > upper_ci | observed < lower_ci, 1, 0)) %>% 
    pull(p) %>% 
    mean()
  
  # Visualize and save plot 
  ggplot(spearman_corr_df, 
         aes(x = dates, 
             y = observed)) +
    geom_line() +
    xlab("Time") +
    ylab("Rank-based Moran's I") +
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
  ggsave(paste0("Figures/RankMoransI/covid_moran_rank_", name, 
                ".pdf", collapse = ""), 
         width = 27, height = 14, unit = "cm")
  
  return(spearman_corr_p)
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
                                     beta_options = list(0, 1, 2, 
                                                         c(1, 1), 
                                                         c(2, 1), 
                                                         c(2, 2), 
                                                         c(1, 1, 1), 
                                                         c(2, 1, 1),
                                                         c(1, 1, 1, 1), 
                                                         c(2, 1, 1, 1),
                                                         c(2, 2, 1, 1),
                                                         c(2, 2, 2, 1), 
                                                         c(1, 1, 1, 1, 1),
                                                         c(2, 1, 1, 1, 1), 
                                                         c(2, 2, 1, 1, 1), 
                                                         c(2, 2, 2, 1, 1)
                                     ), # in form of list
                                     globalalpha = c("TRUE", "FALSE"), 
                                     upper_limit = 5, 
                                     net, 
                                     vts,  
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
                                     
                                     forecast_window) {
  
  train_window <- dim(vts)[1] - forecast_window
  
  # create all possible combinations of alpha and beta order 
  model_options <-  expand.grid(alpha_options, 
                                beta_options, 
                                globalalpha)
  
  model_options$valid <- NA
  
  for(i in seq(1, nrow(model_options))) {
    
    current_option <- model_options[i, ]
    
    model_options[i, ]$valid <-  ifelse(length(current_option$Var2[[1]]) > current_option$Var1 | 
                                          any(current_option$Var2[[1]] > upper_limit),  
                                        FALSE, 
                                        TRUE
                                      )
    
    if(model_options$valid[i]) {
      model_options[i, ]$Var2[[1]] <- c(current_option$Var2[[1]], 
                                        rep(0, current_option$Var1 - length(current_option$Var2[[1]])))
    } else {
      model_options[i, ]$Var2[[1]] <- NA
    }
  }
  
  model_options_valid <- model_options %>% 
    filter(!is.na(Var2), 
           !duplicated(Var2))
  
  # filter out valid parameter combinations 
  # model_options$valid <- (model_options$Var1 == model_options$Var2 %>% 
  #                           lapply(FUN = function(i) {length(i)}) %>% 
  #                           unlist())
  # model_options_valid <- model_options %>% filter(valid == TRUE)
  
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

# Restrictions - tuning parameters ----------------------------------------
# circle through different alpha and beta order combinations and fit 
# GNAR models with the function "fit_and_predict" for each data subset 
fit_and_predict_for_restrictions <- function(alpha_options = seq(1, 5), 
                                             beta_options = list(0, 1, 2, 3, 
                                                                 4, 5, 6, 7, 
                                                                 c(1, 1), 
                                                                 c(2, 1), 
                                                                 c(2, 2), 
                                                                 c(3, 1),
                                                                 c(4, 1),
                                                                 c(5, 1),
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
                                                                 # c(7, 2, 2, 2, 0),
                                                                 c(1, 1, 1, 1, 1),
                                                                 c(2, 1, 1, 1, 1),
                                                                 c(3, 1, 1, 1, 1), 
                                                                 c(4, 1, 1, 1, 1),
                                                                 c(5, 1, 1, 1, 1),
                                                                 c(6, 1, 1, 1, 1),
                                                                 # c(7, 1, 1, 1, 1),
                                                                 c(2, 2, 1, 1, 1), 
                                                                 c(2, 2, 2, 1, 1),
                                                                 c(3, 2, 2, 1, 1),
                                                                 c(4, 2, 2, 1, 1),
                                                                 c(5, 2, 2, 1, 1),
                                                                 c(6, 2, 2, 1, 1),
                                                                 # c(7, 2, 2, 1, 1),
                                                                 c(2, 2, 2, 2, 1), 
                                                                 c(2, 2, 2, 2, 2)
                                             ), # in form of list
                                             globalalpha = TRUE,
                                             upper_limit, 
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
                                             
                                             forecast_window = 5) {
  param <- list()
  
  for (i in seq(1, length(data_list))) {
    res <- fit_and_predict_for_many(vts = data_list[[i]], 
                                    net = net, 
                                    alpha_options = alpha_options, 
                                    beta_options = beta_options, 
                                    globalalpha = globalalpha,
                                    old = old,
                                    forecast_window = forecast_window, 
                                    numeric_vertices = numeric_vertices, 
                                    weight_factor = weight_factor, 
                                    inverse_distance = inverse_distance, 
                                    distance_matrix = distance_matrix, 
                                    weight_index = weight_index,
                                    county_index = county_index, 
                                    upper_limit = upper_limit)
    
    param[[i]] <- return_best_model(results_list = res)
  }
  param_df <- do.call(rbind.data.frame, param)
  
  param_df$data_subset <- seq(1, length(data_list))
  
  return(param_df[, c(3, 1, 2)])
}


save(res, file = "queen_restricted_data_results.RData")




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
  
  param_df <- param_df %>% 
    mutate(order = ifelse(grepl("alpha", 
                                param_df$coefficient_formated), 
                          "alpha-order", 
                          "beta-order"))
  
  
  # visualise the change in coefficient values across data subsets  
  g <- ggplot(param_df, 
                    aes(x = restriction, 
                        y = estimate, 
                        group = coefficient_formated, 
                        color = coefficient_formated)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper, 
                      x = restriction, 
                      group = coefficient_formated, 
                      color = coefficient_formated), 
                  width=0.3, 
                  alpha = 0.5) +
    geom_line(linetype  = "dashed") +
    labs(x = "Pandemic phase", 
         y = "coefficient values", 
         color = "") +
    scale_color_manual(values = c("alpha1" = "#A3A500", 
                                  "alpha2" = "#9590FF", 
                                  "alpha3" = "#FF62BC", 
                                  "alpha4" = "#F8766D", 
                                  "alpha5" = "#E76BF3", 
                                  "beta1.1" = "#00BFC4", 
                                  "beta1.2" = "#00B0F6", 
                                  "beta2.1" = "#00BF7D", 
                                  "beta3.1" = "#39B600", 
                                  "beta4.1" = "#D89000")
                       
    ) +
    theme(legend.position = "bottom") +
    facet_grid(. ~ order, 
               scales = "free")
  
  ggsave(file = paste0("Figures/ParameterDevelopment/", 
                       name,
                       ".pdf", 
                       collapse = ""), 
         plot = g, 
         width = 25, 
         height = 15, unit = "cm")
  
  return(list("plot" = g, 
              "dataset" = param_df))
  
}

# MASE --------------------------------------------------------------------
# compute MASE for certain counties 
compute_MASE <- function(model, 
                         network_name, 
                         n_ahead = 3, 
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
                        counties_subset, 
                        number_counties, 
                        type_name = "grsdt", 
                        types = c("ARIMA", 
                                  "subset_1_gabriel", 
                                  "subset_1_relative", 
                                  "subset_1_soi", 
                                  "subset_1_delaunay", 
                                  "subset_1_train"),
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
                        selected_networks = NULL) {
  # filter current considered network 
  mase_overview_df <- mase_overview %>% 
    filter(type %in% types, 
           CountyName %in% counties_subset) 
  # create factor type 
  mase_overview_df$type <- factor(mase_overview_df$type, 
                                  levels = types)
  
  # filter networks if applicable 
  if(!is.null(selected_networks)) {
    mase_overview_df <- mase_overview_df %>% 
      filter(types %in% selected_networks)
  }
  
  
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
                         type_name = "kdcqeh", 
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
                                         "subset_1_eco_hub" = "#00BF7D"), 
                         selected_networks = NULL) {
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
              selected_networks = selected_networks)
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
                                       type_name = "grsdt", 
                                       counties_subset, 
                                       number_counties, 
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
                                                          "Train")) {
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
                                        type_name = "kdcqeh", 
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
                             type_name = type_name, 
                             counties_subset = counties_subset, 
                             number_counties = number_counties, 
                             color_types = color_types, 
                             label_networks = c("ARIMA", 
                                                "KNN", 
                                                "DNN", 
                                                "Complete",
                                                "Queen", 
                                                "Eco hub"))
}


