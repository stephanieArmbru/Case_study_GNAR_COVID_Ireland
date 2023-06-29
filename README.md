# Instructions 

R code for paper Armbruster and Reinert, 2023, "COVID-19 incidence in the Republic of Ireland:  A case study for network-based time series models"

The structure of the code is modular. <br />
Folder Data includes the following datasets: <br />
• COVID-19_HPSC_County_Statistics_Historic_Data.csv: data set with daily cumulative COVID-19 incidence from Ordnance Survey Ireland, 2023 <br />
• ireland_covid_weekly.csv: data set after pre-processing <br />
• ireland_covid_weekly_2.csv: data with lag-2 COVID-19 incidences after pre-processing <br />
• ireland_shapefile.shp: shapefile for plotting network maps, generated by saving data downloaded from the GADM database <br />
• county_towns_ireland.csv: data set with coordinates for Irish counties from the Ireland Cities Database on simplemaps <br />
• Folder RObjects: contains networks (as GNAR and igraph objects) as well as subsetted data according to COVID-19 regulations <br />

The file functions_paper.R includes all functions written to pre-process data, fit GNAR models, analyse and plot residuals etc. <br />
In order to incorporate the alternative weighting schemes, the functions GNARdesign(), GNARfit() and NofNeighbours() from the package GNAR had to be expanded. The adapted functions GNARdesign_weighting(), GNARfit_weighting() and NofNeighbours_named() are included in the files GNARdesign.R, GNARfit.R and NofNeighbours.R in the folder GNAR. <br />

The file Data_processing.R performs some initial data exploration, the data aggregation to a weekly level and the necessary smoothing for extreme peaks in COVID-19 incidence. <br />
The file Data_processing_restrictions.R processes the data according to COVID-19 restrictions and summarises the datasets belonging to restricted and unrestricted pandemic phases. <br />
The COVID-19 networks are constructed in the file Network_construction.R. <br />
The GNAR models are fit to the restricted and unrestricted datasets, MASE values plotted and residuals testes in the file Model_fitting.R. <br />

The file SM_Model_fitting_entire_datasets.R fits GNAR models to the entire dataset. <br />
The file SM_Model_fitting_subsets.R fits GNAR models to data subsets defined by official COVID-19 guidelines imposed by the Irish Government. <br />

The folder Figures contains all generated figures for the paper. 
