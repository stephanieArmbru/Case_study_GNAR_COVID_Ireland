### PRE-PROCESS COVID-19 DATA FOR IRELAND 

# Install packages --------------------------------------------------------
rm(list=ls())
# set seed to guarantee reproduciblity 
set.seed(1234)

# Load libraries ----------------------------------------------------------
library(readr)
library(magrittr)
library(tidyverse)
library(ISOweek)
library(MASS)
library(xtable) # for tables
library(ade4) # igraph to neighbourhood list object
library(Hmisc) # for weighted variance 
library(Metrics) # for MASE computation 
library(rlist) # for easy concatenation of lists
library(ape)
library(igraph)
library(lubridate) # for dates
library(zoo)

library(latex2exp) # for latex in axis labels 

source("functions_paper.R")
# Load Data ---------------------------------------------------------------
COVID_data_orig <- read_csv("Data/COVID-19_HPSC_County_Statistics_Historic_Data.csv", 
                            show_col_types = FALSE)
COVID_data_orig %>% head()

COVID_data_orig$TimeStamp %>% tail()



# retain the original dataset
COVID_data <- COVID_data_orig %>% 
  as.data.frame() 

# Data exploration --------------------------------------------------------
COVID_data %>% nrow() # 27560 rows
COVID_data %>% ncol() # 16 columns 
COVID_data %>% colnames()

# ObjectID: unique id for each row 
# OrigID: id for county, 26 in total
# CountyName 
# PopulationCensus16: 
# TimeStamp: date and time (separated in preprocessing)
# IGEasting / IGNorthing: geographic identifier
# Lat / Long: gives the "middle" of the county (not its county town)
# UGI: url link
# ConfirmedCovidCases
# PopulationProportionCovidCases
# ConfirmedCovidDeaths
# ConfirmedCovidRecovered
# SHAPE_Length/Area

# check if ObjectID and OrigID identical
COVID_data %>% filter(COVID_data$OBJECTID != COVID_data$ORIGID)
COVID_data %>% dplyr::select(OBJECTID) %>% unique() %>% nrow()
COVID_data %>% dplyr::select(ORIGID) %>% unique() %>% nrow()

# compute relative size of counties  
counties_size <- COVID_data %>% 
  dplyr::select(CountyName, 
                PopulationCensus16) %>% 
  unique() %>% 
  mutate(PopulationChar = as.character(PopulationCensus16), 
         PopulationRel = round(100 * PopulationCensus16 / sum(PopulationCensus16), 
                               digits = 2)) %>% 
  dplyr::select(CountyName, 
                PopulationChar, 
                PopulationRel) %>% 
  arrange(desc(PopulationRel))


# extract time period for which data was recorded  
COVID_data %>% dplyr::select(TimeStamp) %>% head(1) # from 27.02.2020
COVID_data %>% dplyr::select(TimeStamp) %>% tail(1) # to 23.01.2023


# Confirmed COVID-19 Deaths / Recovered
COVID_data %>% 
  dplyr::select(ConfirmedCovidCases) %>% 
  is.na() %>% 
  sum()

COVID_data %>% 
  dplyr::select(ConfirmedCovidDeaths) %>% 
  is.na() %>% 
  sum()
# ConfirmedCovidDeaths is empty column 

COVID_data %>% 
  dplyr::select(ConfirmedCovidRecovered) %>% 
  is.na() %>% 
  sum() 
# ConfirmedCovidRecovered is empty column

COVID_data %>% 
  dplyr::select(PopulationProportionCovidCases) %>% 
  is.na() %>% 
  sum()
# 52 missing values 

COVID_data %>% 
  filter(is.na(PopulationProportionCovidCases)) %>% 
  dplyr::select(TimeStamp) %>% 
  table()
# NA for beginning of pandemic, i.e. between 27.02.2020 and 01.03.2020
# for every county, the data for 27.02 and 01.03 is recorded but NA

COVID_data %>% 
  filter(TimeStamp == "2020/02/28 00:00:00+00")
COVID_data %>% 
  filter(TimeStamp == "2020/02/29 00:00:00+00")
# data for 28.-29.02.2020 is missing  

COVID_data %>% 
  dplyr::select(TimeStamp) %>% 
  unique() %>% 
  head(1)
# first recording on the 27.02.2020

COVID_data %>% 
  filter(PopulationProportionCovidCases != 0, 
         TimeStamp == "2020/03/02 00:00:00+00") %>% 
  dplyr::select(CountyName)
# on 02.03.2020 the first COVID-19 case was recorded in Dublin  


COVID_data %>% 
  group_by(TimeStamp) %>% 
  summarise(count_counties = sum(PopulationProportionCovidCases != 0)) %>% 
  filter(count_counties == 26) %>% 
  head(1)
# on 20.03.2020 COVID-19 was recorded in every county


# Data Pre-processing -----------------------------------------------------
# make object ID and origine ID discrete and 
# remove NA-filled Deaths and Recoveries columns 
COVID_data <- COVID_data %>% 
  mutate(id = as.factor(OBJECTID), 
         origin = as.factor(ORIGID)) %>% 
  dplyr::select(-c(ConfirmedCovidDeaths, ConfirmedCovidRecovered))

# separate time into date and time
COVID_data <- COVID_data %>% 
  separate(col = "TimeStamp", 
           into = c("date", "time"), 
           sep = " ") %>% 
  mutate("date_f" = as.Date(date))

# rank for size 
size_rank <- COVID_data %>%  
  dplyr::select(CountyName, PopulationCensus16) %>% 
  unique() %>%  
  mutate(SizeRank = rank(PopulationCensus16)) %>% 
  dplyr::select(-PopulationCensus16)

COVID_data <- left_join(COVID_data, size_rank, by = "CountyName")

# remove missing dates
COVID_data <- COVID_data %>% 
  drop_na()


# compute daily COVID-19 incidence for 100,000 inhabitants
COVID_data <- COVID_data %>% 
  group_by(CountyName) %>% 
  mutate(dailyCovid = c(0, 
                        diff(PopulationProportionCovidCases, 
                             lag = 1))) %>% 
  ungroup()

# Weekly aggregation ------------------------------------------------------
# aggregate the daily COVID-19 incidence to a weekly COVID-19 incidence 
# for 100,000 inhabitants 
COVID_data_week_agg <- COVID_data %>% 
  dplyr::select(CountyName, 
         date_f, 
         PopulationProportionCovidCases) %>% 
  group_by(CountyName, 
           yw = floor_date(date_f, unit = "week")) %>% 
  summarise(weeklyCasesSum = mean(PopulationProportionCovidCases)) %>% 
  mutate(weeklyCases = c(0, 
                         diff(weeklyCasesSum, 
                              lag = 1))) %>% 
  ungroup()


census <- COVID_data %>% 
  dplyr::select(CountyName, 
         PopulationCensus16) %>% 
  unique()

COVID_data_week_agg <- left_join(x = COVID_data_week_agg, 
                                 y = census, 
                                 by = "CountyName")

COVID_data_week_agg$yw %>% unique %>% length()
# 152 weekly COVID-19 observations for each county 

# plot weekly COVID-19 incidence for each county 
ggplot(COVID_data_week_agg, 
       aes(x = yw, 
           y = weeklyCases, 
           color = CountyName)) + 
  geom_line() +
  xlab("Time") +
  ylab("COVID-19 cases") +
  guides(color = guide_legend(title = "County"))
ggsave("Figures/Visualisation/covid_ireland_weekly_cases.pdf", 
       width = 23, height = 14, unit = "cm")


# Moving average winter 21/22 ----------------------------------------------
# isolate peak in weekly COVID-19 incidence during winter 21/22 and 
# smooth by computing the moving average incidence over 4 consecutive 
# weeks 
COVID_data_week_smooth <- COVID_data_week_agg %>% 
  mutate(weeklyCases = ifelse(yw %in% seq(as.Date("2021-12-12"),
                                          as.Date("2022-02-27"), 
                                          by = 7), # 9 weeks over Christmas
                              rollmean(weeklyCases,
                                        k = 4,
                                        align = 'right',
                                        fill = 0), 
                              weeklyCases))

# plot weekly COVID-19 incidence with smoothed winter season 
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = weeklyCases, 
           color = CountyName)) + 
  geom_line() +
  xlab("Time") +
  ylab("COVID-19-19 cases") +
  guides(color = guide_legend(title = "County"))
ggsave("Figures/Visualisation/covid_ireland_weekly_cases_winter_smoothed.pdf", 
       width = 23, height = 14, unit = "cm")


# Transformation ----------------------------------------------------------
# compute integrated weekly COVID-19 incidence through taking differences
COVID_data_week_smooth <- COVID_data_week_smooth %>%
  group_by(CountyName) %>% 
  mutate(weeklyCasesDiff = c(0, 
                             diff(weeklyCases, 
                                  lag = 1)), 
         weeklyCasesDiff2 = c(0, 
                             diff(weeklyCasesDiff, 
                                  lag = 1))
         ) %>% 
  ungroup()

# identify smallest difference in COVID-19 incidence 
COVID_data_week_smooth[which.min(COVID_data_week_smooth$weeklyCasesDiff), ]$weeklyCasesDiff


# plot difference in weekly COVID-19 incidence 
# stationary and hence used to fit GNAR models 
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = weeklyCasesDiff,
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("1-lag diff. weekly cases") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("Figures/Visualisation/weekly_cases_diff_smoothed.pdf", 
       width = 27, height = 18, unit = "cm")

# plot difference in weekly COVID-19 incidence 
# stationary and hence used to fit GNAR models 
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = weeklyCasesDiff2,
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("1-lag diff. weekly cases") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("Figures/Visualisation/weekly_cases_diff_smoothed_lag_2.pdf", 
       width = 27, height = 18, unit = "cm")


# Box-Cox transformation --------------------------------------------------
# compute Box-Cox transformation for weekly COVID-19 incidence
pdf("Figures/Stationarity/boxcox_weekly_cases.pdf")
boxcox(lm(weeklyCases + 700 ~ 1 + CountyName, 
          data = COVID_data_week_smooth))
# 0 : log transformation 
dev.off()

# compute Box-Cox transformation for 1-lag COVID-19 ID
pdf("Figures/Stationarity/boxcox_weekly_cases_diff.pdf")
boxcox(lm(weeklyCasesDiff + 693 ~ 1 + CountyName, 
          data = COVID_data_week_smooth))
# 1: no transformation 
dev.off()

# Save pre-processed data -------------------------------------------------
# COVID_final <- COVID_data_week_smooth %>% 
#   dplyr::select("CountyName", 
#                 "yw", 
#                 "weeklyCasesSum", 
#                 "weeklyCases", 
#                 "weeklyCasesDiff", 
#                 "PopulationCensus16") %>% 
#   mutate(weeklyCases_non_lag = weeklyCases) %>% 
#   dplyr::select(-weeklyCases) %>% 
#   rename(weeklyCases = weeklyCasesDiff)
# 
# # check if all necessary columns include 
# COVID_final %>% colnames()
# # save data set
# write_csv(COVID_final,
#           file = "Data/ireland_covid_weekly.csv")

