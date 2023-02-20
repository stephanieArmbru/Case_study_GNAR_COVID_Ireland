### SEPARATE DATASET IN SUBSETS ACCORDING TO COVID-19 RESTRICTIONS AND PANDEMIC SITUATIONS

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
COVID_weekly_data <- read_csv(file = "Data/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)

# convert to correct format for GNARfit()
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# convert to data frame with time column for plotting 
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 


# Relevant restrictions ---------------------------------------------------
# create data frame for relevant COVID-19 regulations and their start date
restrictions_df <- data.frame(date = c("28.02.2020", 
                                       "18.08.2020" ,
                                       "26.12.2020" ,
                                       "10.05.2021", 
                                       "06.03.2022", 
                                       "23.01.2023"
) %>% as.Date(format = "%d.%m.%Y"),
content = c("Start - Lockdown", 
            "County-specific lockdowns", 
            "Level-5 lockdown", 
            "Inter-county travel",
            "Ease", 
            "End"
), 
index = seq(1, 6), 
sentiment = c("Restriction", 
              "Restriction", 
              "Restriction",
              "Ease",  
              "Ease",  
              "Ease"
)
)

# create data frame with COVID-19 virus strains and the date at which more 
# than 50% of COVID-19 infections where caused by the strain 
strains <- data.frame(variant = c("Original", 
                                  "Alpha", 
                                  "Delta", 
                                  "Omicron I", 
                                  "Omicron II", 
                                  NA), 
                      onset = c("28/02/2020", 
                                "27/12/2020", 
                                "06/06/2021", 
                                "13/12/2021", 
                                "13/03/2022", 
                                "23/01/2023") %>% 
                        as.Date(format = "%d/%m/%Y"), 
                      index = seq(1, 6))

COVID_weekly_data$restriction <- NA
COVID_weekly_data$restriction_date <- as.Date(NA)
COVID_weekly_data$sentiment <- NA
COVID_weekly_data$strain <- NA

# assign each week the predominant COVID-19 virus strain  
for (i in seq(1, nrow(COVID_weekly_data))) {
  current_date <- COVID_weekly_data$yw[i]
  
  which_strain <- strains[strains$onset > current_date, ]$index %>% min() - 1
  
  COVID_weekly_data$strain[i] <- strains[which_strain, ]$variant
  
  which_restriction <- restrictions_df[restrictions_df$date > current_date, ]$index %>% 
    min() - 1
  
  COVID_weekly_data$phase[i] <- which_restriction
  
  if (which_restriction == 1) {
    which_restriction <- 2
  }
  COVID_weekly_data$restriction[i] <- which_restriction
  COVID_weekly_data$restriction_date[i] <- restrictions_df[which_restriction, ]$date
  COVID_weekly_data$sentiment[i] <- restrictions_df[which_restriction, ]$sentiment
}


# plot restrictions and virus strains 
ggplot(COVID_weekly_data) +
  geom_rect(aes(xmin = yw, xmax = lead(yw),
                ymin = -Inf, ymax = Inf,
                fill = strain)) +
  geom_line(aes(x = yw, 
                y = weeklyCases, 
                group = CountyName, 
                color = phase %>% as.factor())) +
  # geom_vline(aes(xintercept = restriction_date, 
  #                color = restriction %>% as.factor(), 
  #                linetype = sentiment), 
  #            size = 1) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(values = c("Original" = "#FAFAFA", 
                               "Alpha" = "#eeeeee", 
                               "Delta" = "#E1E5E8", 
                               "Omicron I" = "#D0D5D9", 
                               "Omicron II" = "#ABB0B8")) +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave("Figures/Visualisation/covid_id_restrictions.pdf", 
       width = 25, height = 14, unit = "cm")



# Split dataset -----------------------------------------------------------
splits <- c("27.02.2020",  
            "18.08.2020",
            "26.12.2020",
            "10.05.2021",
            "06.03.2022",
            "23.01.2023") %>% 
  as.Date(format = "%d.%m.%Y")

# split COVID-19 data set into separate data sets according to COVID-19 
# regulations 
datasets_list <- list()

for (i in seq(2, length(splits))) {
  datasets_list[[i-1]] <- COVID_weekly_data %>% 
    filter(splits[i-1] < yw & yw <= splits[i]) %>% 
    dplyr::select(CountyName,
                  yw, 
                  weeklyCases) %>% 
    spread(CountyName, 
           weeklyCases) %>% 
    column_to_rownames(var = "yw") %>% 
    as.matrix()
}

# create concatenated data frame with column for pandemic phase 
datasets_plot_list <- list()
for  (i in seq(1, 5)) {
  datasets_plot_list[[i]] <- datasets_list[[i]] %>% 
    as.data.frame() %>% 
    mutate("phase" = i)
  
  datasets_plot_list[[i]]$yw <- rownames(datasets_plot_list[[i]]) %>% as.Date()
}

lapply(datasets_list, FUN = nrow)

datasets_plot <- do.call(rbind.data.frame, datasets_plot_list) %>% 
  gather("county", "weeklyCases", -c(phase, yw))

# plot restriction phases 
ggplot(datasets_plot) +
  geom_line(aes(x = yw, 
                y = weeklyCases, 
                color = phase %>% as.factor(), 
                group = county)) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(breaks = splits %>% as.Date())
ggsave("Figures/Visualisation/covid_id_restrictions_phases.pdf", 
       width = 25, height = 14, unit = "cm")



# Stationarity ------------------------------------------------------------
datasets_plot <- datasets_plot %>% 
  group_by(county) %>% 
  mutate(weeklyCasesDiff2 = c(0, 
                              diff(weeklyCases, 
                                   lag = 1)), 
         weeklyCasesDiff3 = c(0, 
                              diff(weeklyCasesDiff2, 
                                   lag = 1))
  ) %>% 
  ungroup()


# plot restriction phases 
ggplot(datasets_plot) +
  geom_line(aes(x = yw, 
                y = weeklyCasesDiff2, 
                color = phase %>% as.factor(), 
                group = county)) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(breaks = splits %>% as.Date())
ggsave("Figures/Visualisation/covid_id2_restrictions_phases.pdf", 
       width = 25, height = 14, unit = "cm")

# plot restriction phases 
ggplot(datasets_plot) +
  geom_line(aes(x = yw, 
                y = weeklyCasesDiff3, 
                color = phase %>% as.factor(), 
                group = county)) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(breaks = splits %>% as.Date())
ggsave("Figures/Visualisation/covid_id3_restrictions_phases.pdf", 
       width = 25, height = 14, unit = "cm")


# Stationarity for each data set ------------------------------------------
pdf("Figures/Stationarity/boxcox_subset_1.pdf", 
    width = 6, height = 4)
boxcox(lm(COVID_ID ~ CountyName, 
          data = compute_box_cox(1)))
dev.off()

pdf("Figures/Stationarity/boxcox_subset_2.pdf", 
    width = 6, height = 4)
boxcox(lm(COVID_ID ~ CountyName, 
          data = compute_box_cox(2)))
dev.off()

pdf("Figures/Stationarity/boxcox_subset_3.pdf", 
    width = 6, height = 4)
boxcox(lm(COVID_ID ~ CountyName, 
          data = compute_box_cox(3)))
dev.off()

pdf("Figures/Stationarity/boxcox_subset_4.pdf", 
    width = 6, height = 4)
boxcox(lm(COVID_ID ~ CountyName, 
          data = compute_box_cox(4)))
dev.off()

pdf("Figures/Stationarity/boxcox_subset_5.pdf", 
    width = 6, height = 4)
boxcox(lm(COVID_ID ~ CountyName, 
          data = compute_box_cox(5)))
dev.off()


# Variability for each data set --------------------------------------------
# compute variance in 1-lag COVID-19 ID for each data subset
lapply(datasets_list, FUN = function(i) {
  apply(i, MARGIN = 2, FUN = function(j) {sd(j)}) %>% 
    mean() %>% 
    round(2)
})

# number of weeks in every data subset
lapply(datasets_list, FUN = function(i) {
  nrow(i)
})



# Pandemic situation ------------------------------------------------------
data_restrictive <- rbind(datasets_list[[1]], 
                          datasets_list[[3]]) %>% 
  as.data.frame()
data_restrictive %>% nrow()

# add missing data
data_restrictive_NA <- data_restrictive %>% 
  mutate(time = as.Date(rownames(data_restrictive))) %>% 
  complete(time = seq.Date(min(time), max(time), by="week")) 

times_restrictive <- data_restrictive_NA$time
data_restrictive_NA <- data_restrictive_NA %>% 
  dplyr::select(-time)

rownames(data_restrictive_NA) <- seq.Date(min(times_restrictive), 
                                          max(times_restrictive), 
                                          by="week")

# which dates are missing 
setdiff(rownames(data_restrictive_NA), rownames(data_restrictive))

# predicted dates
rownames(data_restrictive) %>% tail(5)



data_free <- rbind(datasets_list[[2]], 
                   datasets_list[[4]],
                   datasets_list[[5]]) %>% 
  as.data.frame()
data_free %>% nrow()

# add missing data 
data_free_NA <- data_free %>% 
  mutate(time = as.Date(rownames(data_free))) %>% 
  complete(time = seq.Date(min(time), max(time), by="week"))
times_free <- data_free_NA$time
data_free_NA <- data_free_NA %>% 
  dplyr::select(-time)

rownames(data_free_NA) <- seq.Date(min(times_free), 
                                   max(times_free), 
                                   by="week")
# which dates are missing
setdiff(rownames(data_free_NA), rownames(data_free))

# predicted dates
rownames(data_free) %>% tail(5)



datasets_list_coarse <- list("restrictive" = data_restrictive_NA %>% 
                               as.matrix(), 
                             "free" = data_free_NA %>% 
                               as.matrix())


data_plot <- rbind(data_restrictive %>% 
                     mutate(phase = "restrictive", 
                            yw = rownames(data_restrictive) %>% 
                              as.Date()), 
                   data_free %>% mutate(phase = "free", 
                                        yw = rownames(data_free) %>% 
                                          as.Date())) %>% 
  gather("county", "weeklyCases", -c(phase, yw)) 


# plot restriction phases 
ggplot(data_plot) +
  geom_line(aes(x = yw, 
                y = weeklyCases, 
                color = phase, 
                group = county)) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(breaks = splits %>% as.Date())
ggsave("Figures/Visualisation/covid_id_pandemic_phases.pdf", 
       width = 25, height = 14, unit = "cm")


# standard deviation for subsets 
lapply(datasets_list_coarse, FUN = function(i) {
  apply(i, MARGIN = 2, FUN = function(j) {sd(j %>% na.omit())}) %>% 
    mean() %>% 
    round(2)
})

# Save data subsets -------------------------------------------------------
save(datasets_list,
     file = "Data/RObjects/data_subsets_regulations.RData")
save(datasets_list_coarse, 
     file = "Data/RObjects/data_subsets_pandemic_situations.RData")


