## Load packages

library(LakeEnsemblR)
library(rLakeAnalyzer)
library(ggplot2)
library(patchwork)
library(dplyr)
library(lubridate)


## load data on lakes
adkl <- read.csv("adk_lakes4sim.csv")
## load meteorological driver data
era5 <- arrow::read_parquet("era5_adk_1980-2022.parquet") %>% data.table::as.data.table()
## create list of grid cells
gridcells <- expand.grid(long = unique(era5$lon), lat = unique(era5$lat))


## compute hypsographic curve information for each lake
allbathy <- lapply(1:nrow(adkl), function(x){
  approx.bathy(ceiling(adkl$mdep[x]), adkl$area_ha[x] * 10000, adkl$meandep[x], method = "voldev")
})


# Run simulation for each lake

## set yaml file identity
ler_yaml <- "LakeEnsemblR.yaml"

t0 <- Sys.time() # save start time
for(i in 1:nrow(adkl)){  # initiate loop through lakes
  
  # get location for lake i
  dm1 <- c(long = adkl$lon[i], lat = adkl$lat[i])
  
  # identify which grid cell is closest to lake
  era5gridno <- gridcells[which.min(as.matrix(dist(rbind(dm1, gridcells)))[-1,1]),]
  
  # select weather from appropriate grid
  era5 %>% filter(lat == era5gridno$lat, lon == era5gridno$long) %>% 
    select(-lat, -lon) %>% 
    # filter(year(datetime) %in% 1980:2022) %>% 
    arrange(datetime) %>% 
    mutate(datetime = as.character(format(datetime))) %>% 
    write.csv("ler_adk_met.csv", row.names = FALSE) # save weather driver to common filename
  
  # set filename for weather driver data in yaml
  input_yaml_multiple(file = ler_yaml, "ler_adk_met.csv", key1 = "input", key2 = "meteo", key3 = "file")
  
  # set parameters in yaml file
  ## bathymetry
  hyps <- allbathy[[i]]
  colnames(hyps) <- c("Depth_meter", "Area_meterSquared")
  write.csv(hyps, 'LakeEnsemblR_bathymetry_standard.csv', row.names = FALSE)
  ## initial temps of 5 C
  write.csv(data.frame(depths = allbathy[[i]]$depths, temps = 5), "init_temps.csv", row.names = FALSE)
  
  ## identify max depth
  maxdep <- max(allbathy[[i]]$depths)
  
  ## set max and intial depth to max depth in yaml file
  input_yaml_multiple(file = ler_yaml, maxdep, key1 = "location", key2 = "depth")
  input_yaml_multiple(file = ler_yaml, maxdep, key1 = "location", key2 = "init_depth")
  
  ## set kd
  ### https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2427.2000.00518.x
  # input_yaml_multiple(file = ler_yaml, 0.15 * adkl$meanvalue_doc[i]^1.08,
  #                     key1 = "input", key2 = "light", key3 = "Kw")
  
  ## compute kw based on trend in secchi depth
  kwtrend <- 1.7/(adkl$sec[i] + (adkl$trend[i] * 0:42))
  ## force kw to be between 0 and 5.5 to remove errors caused by linear trends
  kwtrend[kwtrend > 5.5 | kwtrend < 0] <- 5.66
  
  ## save time varying kw values
  write.csv(data.frame(datetime = seq(ymd_hms("1980-01-01 12:00:00"),ymd_hms("2022-12-31 12:00:00"), "1 year"), 
                       Extinction_coefficient_perMeter = kwtrend), 
            "kdvals.csv", row.names = FALSE)
  ## set time varying kw value filename to yaml file
  input_yaml_multiple(file = ler_yaml, "kdvals.csv",
                      key1 = "input", key2 = "light", key3 = "Kw")
  
  
  ## Assume lake is a circle
  ### set fetch, length, and width parameters as the diameter
  fetch <- sqrt((adkl$area_ha[i] * 10000)/pi)*2 
  blen <- sqrt((adkl$area_ha[i] * 10000)/pi)*2
  bwid <- sqrt((adkl$area_ha[i] * 10000)/pi)*2 
  
  ### set fetch, length, and width parameters in yaml file
  input_yaml_multiple(file = ler_yaml, fetch, key1 = "model_parameters", key2 = "FLake", key3 = "fetch_lk")
  input_yaml_multiple(file = ler_yaml, blen, key1 = "model_parameters", key2 = "GLM", key3 = "bsn_len")
  input_yaml_multiple(file = ler_yaml, bwid, key1 = "model_parameters", key2 = "GLM", key3 = "bsn_wid")
  
  ## set elevation in yaml file
  input_yaml_multiple(file = ler_yaml, adkl$elevation_est[i], key1 = "location", key2 = "elevation")
  
  # set start and end times in yaml file
  input_yaml_multiple(file = ler_yaml, "1980-01-01 01:00:00", key1 = "time", key2 = "start")
  input_yaml_multiple(file = ler_yaml, "2022-12-31 23:00:00", key1 = "time", key2 = "stop")
  
  # set output file id
  input_yaml_multiple(file = ler_yaml, paste0("adk_trend/output_", adkl$Permanent_[i]), key1 = "output", key2 = "file")
  
  # setup simulation config
  config_file <- ler_yaml
  model <- c("GLM", "GOTM", "Simstrat")
  
  export_config(config_file = config_file, model = model)
  
  # run model ensemble
  run_ensemble(config_file = config_file, model = model, parallel = FALSE)
  
}
t1 <- Sys.time()