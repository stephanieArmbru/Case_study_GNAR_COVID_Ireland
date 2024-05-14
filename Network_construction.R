### NETWORK CONSTRUCTION

# Install packages --------------------------------------------------------
rm(list=ls())
# set seed to guarantee reproducibility 
set.seed(1234)

# Load libraries ----------------------------------------------------------
library(ade4) # igraph to neighbourhood list object
library(Hmisc) # for weighted variance 
library(Metrics) # for MASE computation 
library(rlist) # for easy concatenation of lists
library(ape)
library(readr)
library(igraph)
library(GNAR)
library(tidyverse)
library(xtable) # for tables 
library(geosphere) # for long / lat distance calculation

# for MAPS: 
library(sp)
library(raster) 
library(leaflet)
library(mapview) # for saving

# for neighbourhood construction 
library(spdep) 
# to create neighbourhood data frame
# library(expp) 

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


# Load functions 
source(file = "functions_paper.R")

# Turn off warnings
options(warn = -1)

# Set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load pre-processed data -------------------------------------------------
# weekly data
COVID_weekly_data <- read_csv(file = "Data/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)
# data for pandemic phases
load(file = "Data/RObjects/data_subsets_pandemic_situations.RData")

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


# Parameters --------------------------------------------------------------
# Economic hubs
hubs <- c("Dublin", 
          "Cork",
          "Limerick", 
          "Galway", 
          "Waterford")

# Load maps data ----------------------------------------------------------
# download Ireland data level 1 from GADM database
ireland <- raster::getData('GADM',
                           country='IRL',
                           level = 1)

# comparison of county spelling in maps data and COVID data 
setdiff(COVID_weekly_data$CountyName %>% unique(), ireland$NAME_1)
setdiff(ireland$NAME_1, COVID_weekly_data$CountyName %>% unique())
# different writing: Laois and Laoighis (the former is the newer spelling)

# correct spelling in ireland data
ireland$NAME_1[ireland$NAME_1 == "Laoighis"] <- "Laois"

# save Spatial Polygon Data Frame as shapefile 
shapefile(x = ireland,
          file = "Data/ireland_shapefile.shp", 
          overwrite = TRUE)


# read Ireland shapefile 
ireland_shp <- st_read("Data/ireland_shapefile.shp")

# construct centroid for each county 
coord <- ireland_shp %>%
  st_geometry() %>%
  st_centroid() %>%
  st_coordinates()

# assign county names to their centroid coordinates 
rownames(coord) <- ireland_shp$NAME_1


# read data frame including county towns and their coordinates 
county_towns <- read_delim("Data/county_towns_ireland.csv", 
                           delim = ";", 
                           escape_double = FALSE, 
                           trim_ws = TRUE, 
                           show_col_types = FALSE, 
                           locale = locale(decimal_mark = ","))

# match county towns and COVID incidence 
c1 <- county_towns %>% dplyr::pull(admin_name) %>% unique()
c2 <- COVID_weekly_data$CountyName %>% unique()

# check if all counties are included 
diff_c <- setdiff(c1, c2)
same_c <- intersect(c1, c2)

# filter out all relevant county towns 
county_towns_lim <- county_towns %>% 
  filter(admin_name %in% same_c)

# investigate counties with multiple county towns 
double_county <- county_towns_lim %>% 
  dplyr::pull(admin_name) %>% 
  duplicated()

county_towns_lim %>% 
  filter(admin_name %in% 
           county_towns_lim[double_county, ]$admin_name)

# manually set county towns for counties with multiple county towns 
county_towns_lim_unique <- county_towns_lim %>% 
  filter(city != "Drogheda", 
         city != "Nenagh")

# check if right number of counties included 
county_towns_lim_unique %>% nrow()

# construct matrix with coordinates 
coord_central_towns <- county_towns_lim_unique %>% 
  dplyr::select(lat, lng)

# copy data frame for two different rownames 
coord_central_towns_counties <- coord_central_towns

# assign county towns as rownames 
rownames(coord_central_towns) <- county_towns_lim_unique %>% 
  pull(city)
# assign county names as rownames 
rownames(coord_central_towns_counties) <- county_towns_lim_unique %>% 
  pull(admin_name)

# transform into matrix
coords_ct <- as.matrix(coord_central_towns)[, c(2,1)]
coords_ct_counties <- as.matrix(coord_central_towns_counties)[, c(2,1)]


# Classification ----------------------------------------------------------
# classification of counties into urban and rural according to Census 2016 data
urban <- c("Dublin", "Louth", "Meath", 
           "Kildare", "Wicklow", "Waterford", 
           "Limerick", "Galway", "Cork", "Longford")
rural <- c("Wexford", "Kilkenny", "Carlow", 
           "Laois", "Kerry", "Tipperary", 
           "Offaly", "Westmeath", "Cavan", 
           "Monaghan", "Leitrim", "Mayo", 
           "Sligo", "Donegal", "Clare",
           "Roscommon")

# check that no counties are double classified, and all counties are included 
(rural %in% urban) %>% any()
c(rural, urban) %>% length()

colnames(covid_cases) %in% c(rural, urban) %>% all()

# create county name vector in the correct order  
counties <- covid_cases %>% colnames()

# create factor vector indicating urbanisation status for each county
urbanisation_factor <- c()
for (name in counties) {
  if (name %in% rural) {
    urbanisation_factor <- c(urbanisation_factor, "rural")
  } 
  if (name %in% urban) {
    urbanisation_factor <- c(urbanisation_factor, "urban")
  }
}

# reformat as factor 
urbanisation_factor <- urbanisation_factor %>% as.factor()


# create mixture coordinate matrix 
# rural: centroid coordinates as node representation
# urbal: county town coordinates as node representation 

# guarantee same ordering 
coord_ord <- coord[match(counties, rownames(coord)), ]
coords_ct_counties_ord <- coords_ct_counties[match(counties, 
                                                   rownames(coords_ct_counties)), 
]

# filter coordinates for each county
coord_rural <- coord_ord[urbanisation_factor == "rural", ]
coord_urban <- coords_ct_counties_ord[urbanisation_factor == "urban", ]

# create matrix with coordinates for every county 
coord_urbanisation_prelim <- rbind(coord_rural, 
                                   coord_urban)
# guarantee same ordering 
coord_urbanisation <- coord_urbanisation_prelim[match(counties, 
                                                      rownames(coord_urbanisation_prelim)), 
]

# check if order is correct 
all(rownames(coord_urbanisation) == counties)

# Great Circle distance ---------------------------------------------------
# compute pairwise Great Circle distances between all counties 
dist_urbanisation <- circle_distance(coord_urbanisation)

# construct data frame including population density in Tsd. for each county
population_weight <- COVID_weekly_data %>% 
  dplyr::select(CountyName, 
                PopulationCensus16) %>% 
  unique() %>% 
  mutate(weight = PopulationCensus16 / 1000) %>% 
  dplyr::select(-PopulationCensus16)

# General maps ------------------------------------------------------------
# plot maps of counties
pal <- colorFactor("Reds", ireland$NAME_1)

leaflet(ireland) %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(NAME_1),
              highlightOptions = highlightOptions(color = "white", weight =2,
                                                  bringToFront = TRUE),
              label = ~NAME_1)

# visualise centres according to urbanisation state 
# plot(st_geometry(ireland_shp),
#      border="grey",
#      col = c("#78BE21", "#516FC4")[urbanisation_factor],
#      legend = c("rural", "urban"))
# points(coord_urbanisation[, 1],
#        coord_urbanisation[, 2])
# addTextLabels(xCoords = coord_urbanisation[, 1],
#               yCoords = coord_urbanisation[, 2],
#               labels=rownames(coord_urbanisation),
#               cex.label = 0.9,
#               col.label = "black",
#               col.line = "black")


# Network construction ----------------------------------------------------
# 1. Queen's contiguity: connection to all neighbouring counties 
nb_list_queen <- poly2nb(ireland, 
                         queen = TRUE, 
                         row.names = ireland$NAME_1)

# extract neighbourhood (nb) list and convert to igraph object 
covid_net_queen_igraph <- neighborsDataFrame(nb = nb_list_queen) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() # simplify to obtain no self-loops and no multiple edges

# create GNAR object 
covid_net_queen_gnar <- covid_net_queen_igraph %>% igraphtoGNAR()

# create ordered county index data frame
county_index_queen <- data.frame("CountyName" = covid_net_queen_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(nb_list_queen,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add = TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8,
#      font = 2,
#      pos = 1)

# 2. Rook's contiguity
nb_list_rook <- poly2nb(ireland, 
                        queen = FALSE, 
                        row.names = ireland$NAME_1)

# analyse differences between Queen's and Rook's network
diffnb(x = nb_list_queen, 
       y = nb_list_rook)
# produces same links due to non-grid like structure 

# 3. Economic hubs
# network with edges to the direct neighbors and nearest economic hub
covid_net_eco_hubs <- create_eco_hub_igraph(dist_urbanisation, 
                                            coord_urbanisation)
# create igraph object
covid_net_eco_hubs_igraph <- covid_net_eco_hubs$igraph_net
# create GNAR object
covid_net_eco_hubs_gnar <- covid_net_eco_hubs_igraph %>% igraphtoGNAR()

# extract nb list from igraph object 
nb_eco_hubs <- igraph2nb(gr = covid_net_eco_hubs_igraph)

ord_coord_eco_hubs <- covid_net_eco_hubs$ordered_coord
coord_hubs <- covid_net_eco_hubs$ordered_hubs

# create county index data frame 
county_index_eco_hubs <- data.frame("CountyName" = ord_coord_eco_hubs %>%
                                      rownames(), 
                                    "index" = seq(1, 26))

# visualise network 
# plot(st_geometry(ireland_shp),
#      border = "grey")
# plot(nb_eco_hubs,
#      ord_coord_eco_hubs,
#      add = TRUE,
#      pch = 19, cex = 0.6)
# text(coord_hubs$X,
#      coord_hubs$Y,
#      labels = rownames(coord_hubs),
#      cex = 0.8, font = 2, pos = 1,
#      col = "#516FC4")


# 3. Railway-based network 
# Newbridge is biggest town in Kildare, Naas is county town, 
train_connections <- matrix(c("Dublin", "Dún Dealgan", 
                              "Dublin", "Trim", 
                              "Dublin", "Naas", 
                              "Dublin", "Wicklow", 
                              "Cork", "Tralee", 
                              "Cork", "Limerick", 
                              "Cork", "Clonmel", 
                              "Galway", "Ennis", 
                              "Galway", "Tullamore", 
                              "Galway", "Ros Comáin", 
                              "Limerick", "Ennis", 
                              "Limerick", "Clonmel", 
                              "Limerick", "Port Laoise", 
                              "Limerick",   "Cork", 
                              "Waterford", "Kilkenny", 
                              "Waterford", "Clonmel", 
                              "Dún Dealgan", "Dublin", # Dundalk 
                              "Tralee", "Cork", 
                              "Carlow", "Naas", 
                              "Carlow", "Kilkenny", 
                              "Carlow", "Tullamore", 
                              "Carlow", "Port Laoise",
                              "Ennis", "Galway", 
                              "Ennis",    "Limerick", 
                              "Kilkenny", "Carlow", 
                              "Kilkenny", "Waterford", 
                              "Naas", "Dublin", 
                              "Naas", "Carlow", 
                              "Naas",   "Port Laoise", 
                              "Naas",   "Tullamore", # Newbridge
                              "Sligo", c("Carrick on Shannon"), 
                              "Ros Comáin", "Castlebar", 
                              "Ros Comáin",  "Tullamore", 
                              "Ros Comáin", "Galway", # Roscommon
                              "Mullingar", "Trim", 
                              "Mullingar", "Longford", 
                              "Wicklow", "Dublin", 
                              "Wicklow", "Wexford", 
                              "Clonmel", "Waterford", 
                              "Clonmel", "Limerick", 
                              "Clonmel", "Cork", 
                              "Clonmel", "Tralee", 
                              "Clonmel", "Port Laoise",
                              "Wexford", "Wicklow", 
                              "Longford", "Carrick on Shannon", 
                              "Longford", "Mullingar", 
                              "Trim", "Mullingar", 
                              "Trim", "Dublin", # Enfield
                              "Carrick on Shannon", "Sligo", 
                              "Carrick on Shannon", "Longford", 
                              "Tullamore", "Ros Comáin", 
                              "Tullamore", "Port Laoise", 
                              "Tullamore", "Naas", 
                              "Tullamore", "Carlow", 
                              "Tullamore", "Galway", 
                              "Port Laoise", "Tullamore", 
                              "Port Laoise", "Naas",
                              "Port Laoise", "Carlow", 
                              "Port Laoise", "Limerick", 
                              "Port Laoise", "Cork", 
                              "Port Laoise", "Tralee", 
                              "Castlebar", "Ros Comáin", 
                              "Monaghan", "Dún Dealgan", 
                              "An Cabhán", "Longford", 
                              "Lifford", "Sligo"
), ncol = 2, byrow = TRUE) %>% 
  as.data.frame()
colnames(train_connections) <- c("city", "city2")


no_train_counties <- c("Cavan", 
                       "Donegal", 
                       "Monaghan")
# "An Cabhán" = NA, # Cavan, not reachable by train 
# "Monaghan" = NA, # not reachable by train 
# "Lifford" = NA # Donegal, not reachable by train


# substitute county towns with county index for nb list generation
county_to_town <- county_towns_lim_unique %>% 
  dplyr::select("city", "admin_name") %>%
  mutate("name_index" = seq(1, 26))

# save order to order coordinates adequately for plot 
order_names <- county_to_town %>% 
  arrange(name_index) %>% 
  pull(admin_name)

town_to_index <- county_to_town %>% 
  dplyr::select(-"admin_name")

# substitute town names with corresponding indices since function edgelist() 
# requires numeric vertices
first_step <- left_join(train_connections, town_to_index, by = "city")
colnames(first_step) <- c("sub", "city", "county1")
second_step <- left_join(first_step, town_to_index, by = "city")
colnames(second_step) <- c("sub", "sub2", "county1", "county2")

# create edge list
train_connections_m <- second_step %>% 
  dplyr::select("county1",
                "county2") %>% 
  as.matrix(ncol = 2, byrow = TRUE)

# create igraph object from edge list
covid_net_train_igraph <- graph_from_edgelist(train_connections_m, 
                                              directed = FALSE)
# create GNAR object 
covid_net_train_gnar <- covid_net_train_igraph %>% igraphtoGNAR()

# create nb list object 
nb_list_train <- igraph2nb(gr = covid_net_train_igraph)

# reorder so the ordering of the nb list and the coordinates match
coord_urb_ord <- coord_urbanisation[match(order_names, 
                                          rownames(coord_urbanisation)), ] %>% 
  as.data.frame()

# create ordered county index data frame 
county_index_train <- data.frame("CountyName" = coord_urb_ord %>% rownames(), 
                                 "index" = seq(1, 26))

# visualise network 
# plot(st_geometry(ireland_shp),
#      border = "grey")
# plot(nb_list_train,
#      coord_urb_ord,
#      add = TRUE,
#      pch = 19, cex = 0.6)
# text(coord_urb_ord$X,
#      coord_urb_ord$Y,
#      labels = rownames(coord_urb_ord),
#      cex = 0.8, font = 2, pos = 1)


# 4. Delaunay triangulation 
nb_list_delaunay <- tri2nb(coords = coord_urbanisation, 
                           row.names = coord_urbanisation %>% 
                             rownames())

# create igraph object
covid_net_delaunay_igraph <- neighborsDataFrame(nb = nb_list_delaunay) %>% 
  graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify()

# create GNAR object 
covid_net_delaunay_gnar <- covid_net_delaunay_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_delaunay <- data.frame("CountyName" = covid_net_delaunay_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))
# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(nb_list_delaunay,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)


# 5. Gabriel neighbourhoods 
nb_list_gabriel <- gabrielneigh(coords = coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames(), sym = TRUE)
# symmetry required for every region to have edges 

# create igraph object
covid_net_gabriel_igraph <- neighborsDataFrame(nb = nb_list_gabriel) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object
covid_net_gabriel_gnar <- covid_net_gabriel_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_gabriel <- data.frame("CountyName" = covid_net_gabriel_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# visualise network
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(nb_list_gabriel,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

# 6. Relative neighbourhood
nb_list_relative <- relativeneigh(coords = coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames(), 
           sym = TRUE)

# create igraph object 
covid_net_relative_igraph <- neighborsDataFrame(nb = nb_list_relative) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
covid_net_relative_gnar <- covid_net_relative_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_relative <- data.frame("CountyName" = covid_net_relative_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(nb_list_relative,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

# 7. Sphere of Influence (SOI) neighbourhood 
nb_list_soi <- soi.graph(nb_list_delaunay, 
                         coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames())

# create igraph object 
covid_net_soi_igraph <- neighborsDataFrame(nb = nb_list_soi) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
covid_net_soi_gnar <- covid_net_soi_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_soi <- data.frame("CountyName" = covid_net_soi_igraph %>%
                                 V() %>% 
                                 names(), 
                               "index" = seq(1, 26))

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(nb_list_soi,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)


# generate Complete network as KNN network with the maximum number of neighbours
complete_net <- knearneigh(x = coord_urbanisation, 
                           k = 25, 
                           longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames())

# create igraph object 
complete_net_igraph <- neighborsDataFrame(complete_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
complete_net_gnar <- complete_net_igraph %>% 
  igraphtoGNAR()

# create order county index data frame 
county_index_complete <- data.frame("CountyName" = complete_net_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))


# KNN and DNN networks (best for restricted and unrestricted data sets)
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

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(knn_11,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

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

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(knn_21,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

# DNN (d = 125)
dnn_125 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 200,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_125_igraph<- neighborsDataFrame(nb = dnn_125) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_125_gnar <- dnn_125_igraph %>% 
  igraphtoGNAR()

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(dnn_125,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)

# DNN (d = 325)
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

# visualise network 
# plot(st_geometry(ireland_shp),
#      border="grey")
# plot(dnn_325,
#      coord_urbanisation,
#      pch = 19, cex = 0.6,
#      add=TRUE)
# text(coord_urbanisation[, 1],
#      coord_urbanisation[, 2],
#      labels = rownames(coord_urbanisation),
#      cex = 0.8, font = 2, pos = 1)


# Save objects -----------------------------------------------------------
# save neighbourhood lists 
save(list = "nb_list_queen", 
     file = "Data/RObjects/nb_list.RData")

# save GNAR objects for every network
save(list = c("covid_net_queen_gnar",
              "covid_net_eco_hubs_gnar",
              "covid_net_train_gnar",
              "covid_net_delaunay_gnar",
              "covid_net_gabriel_gnar",
              "covid_net_relative_gnar",
              "covid_net_soi_gnar",
              "complete_net_gnar"),
     file = "Data/RObjects/GNAR.RData")

# save igraph objects for every network
save(list = c("covid_net_queen_igraph",
              "covid_net_eco_hubs_igraph",
              "covid_net_train_igraph",
              "covid_net_delaunay_igraph",
              "covid_net_gabriel_igraph",
              "covid_net_relative_igraph",
              "covid_net_soi_igraph",
              "complete_net_igraph"),
     file = "Data/RObjects/igraph.RData")

# save county indices for every network
save(list = c("county_index_queen",
              "county_index_eco_hubs",
              "county_index_train",
              "county_index_delaunay",
              "county_index_gabriel",
              "county_index_relative",
              "county_index_soi",
              "county_index_complete"),
     file = "Data/RObjects/county_index.RData")

# save vector with population sizes
save(population_weight,
     file = "Data/RObjects/population_weight.RData")

# save matrix of distances between counties
save(dist_urbanisation,
     file = "Data/RObjects/distance_urbanisation.RData")

# save coordinates to represent counties 
# (centroids for rural, county towns for urban)
save(coord_urbanisation,
     file = "Data/RObjects/coord_urbanisation.RData")

save(urbanisation_factor, 
     file = "Data/RObjects/urbanisation_factor.RData")

