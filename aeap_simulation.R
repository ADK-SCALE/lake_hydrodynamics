library(LakeEnsemblR)
library(rLakeAnalyzer)
library(ggplot2)
library(patchwork)
library(adklakedata)
library(dplyr)


t1 <- Sys.time()
# set file name for ler yaml
ler_yaml <- "LakeEnsemblR.yaml"

# only look at lakes with max depth greater than 5 m
meta <- adk_data("meta")
meta <- meta[meta$max.depth > 5,]

era5 <- arrow::read_parquet("era5_adk_1980-2022.parquet") %>% 
  data.table::as.data.table()

tdo <- adk_data("tempdo")

meta2 <- read.csv("C:/Users/borre/Desktop/JP/adk.csv")
meta2 <- meta2[meta2$max.depth > 5,]


# get secchi depths
secchi <- adk_data("secchi")

t1 <- Sys.time()

for(i in 1:nrow(meta)){
  
  rn <- which.min(as.matrix(dist(rbind(data.frame(lat = meta$lat[i], lon = meta$long[i]),
                                       unique(dplyr::select(era5, lat, lon)))))[-1,1])
  
  latlon <- unique(dplyr::select(era5, lat, lon))[rn,]
  
  ## Set location data (from `meta`)
  input_yaml_multiple(file = ler_yaml, meta$lake.name[i], key1 = "location", key2 = "name")
  input_yaml_multiple(file = ler_yaml, meta$lat[i], key1 = "location", key2 = "latitude")
  input_yaml_multiple(file = ler_yaml, meta$long[i], key1 = "location", key2 = "longitude")
  input_yaml_multiple(file = ler_yaml, meta$elevation.m[i], key1 = "location", key2 = "elevation")
  input_yaml_multiple(file = ler_yaml, meta$max.depth[i], key1 = "location", key2 = "depth")
  
  
  ## Set hypsographic curve by approximating bathymetry using rLakeAnalyzer
  bathy <- approx.bathy(meta$max.depth[i], meta$SA.ha[i]*1e5, meta$mean.depth[i])
  colnames(bathy) <- c("Depth_meter", "Area_meterSquared")
  write.csv(bathy, "lake_hyps.csv", row.names = FALSE)
  
  ## set hypsography and max depth in yaml
  input_yaml_multiple(file = ler_yaml, "lake_hyps.csv", key1 = "location", key2 = "hypsograph")
  input_yaml_multiple(file = ler_yaml, meta$max.depth[i], key1 = "location", key2 = "init_depth")
  
  
  # compute initial profile
  laketdo <- filter(tdo, lake.name == meta$lake.name[i])
  initlake <- filter(laketdo, date == min(laketdo$date), depth <= meta$max.depth[i]) %>% select(deps = depth, temp)
  initlake$temp <- 4
  write.csv(initlake, "init_profile.csv", row.names = FALSE)
  
  ## set initial prof in yaml
  input_yaml_multiple(file = ler_yaml, "init_profile.csv", key1 = "input", key2 = "init_temp_profile", key3 = "file")
  
  ## compute light ext
  kdval <- secchi %>% group_by(lake.name) %>% 
    summarize(kd = median(1.7/secchi, na.rm = TRUE))
  
  kdval2 <- secchi %>% mutate(kd = 1.7/secchi) %>% 
    filter(lake.name == meta$lake.name[i]) %>% 
    select(date, kd) %>% mutate(kd = zoo::na.aggregate(kd), kd= round(kd, 3)) %>% 
    arrange(date) %>% 
    rename(datetime = date, Extinction_Coefficient_perMeter = kd) %>% 
    mutate(datetime = paste(datetime, "12:00:00"))
  kdval2 <- rbind(kdval2, data.frame(datetime = c("1980-01-01 12:00:00", "2022-12-31 12:00:00"), 
                                     Extinction_Coefficient_perMeter = round(kdval$kd[kdval$lake.name %in% meta$lake.name[i]], 3))) %>%
    arrange(datetime)
  
  
  write.csv(kdval2, "kdvals.csv", row.names = FALSE)
  ## set kd in yaml
  # input_yaml_multiple(file = ler_yaml, round(kdval$kd[kdval$lake.name %in% meta$lake.name[i]], 2),
  #                     key1 = "input", key2 = "light", key3 = "Kw")
  input_yaml_multiple(file = ler_yaml, "kdvals.csv",
                      key1 = "input", key2 = "light", key3 = "Kw")
  
  ## set start date
  input_yaml_multiple(file = ler_yaml, "1980-01-01 00:00:00", key1 = "time", key2 = "start")
  
  ## set end date
  input_yaml_multiple(file = ler_yaml, "2022-12-31 00:00:00", key1 = "time", key2 = "stop")
  
  ## Set observations
  ler_obs <- tdo %>% filter(lake.name == meta$lake.name[i], !is.na(temp)) %>% 
    select(datetime = date, Depth_meter = depth, Water_Temperature_celsius = temp)
  
  
  write.csv(ler_obs, "LakeEnsemblR_wtemp_profile_standard.csv", row.names = FALSE)
  
  input_yaml_multiple(file = ler_yaml, "LakeEnsemblR_wtemp_profile_standard.csv", key1 = "observations", key2 = "temperature", key3 = "file")
  
  
  ## Set meteorology driver file
  lakemets <- era5 %>% filter(lat == latlon$lat[1], lon == latlon$lon[1]) %>% select(-lat, - lon) %>% arrange(datetime)
  
  lakemets$datetime <- as.character(format(lakemets$datetime))
  write.csv(lakemets, "LakeEnsemblr_meteo_adk_sub.csv", row.names = FALSE)
  
  # input_yaml_multiple(file = ler_yaml, "LakeEnsemblr_meteo_narr_sub.csv", key1 = "meteo", key2 = "file")
  input_yaml_multiple(file = ler_yaml, "LakeEnsemblr_meteo_adk_sub.csv", key1 = "meteo", key2 = "file")
  
  ## Set output file
  outfile <- paste0("output_ice_", gsub(" ", "", meta$lake.name[i]))
  input_yaml_multiple(file = ler_yaml, outfile, key1 = "output", key2 = "file")
  
  # params
  ## use length/width measurements
  ## Assume circle (commented out)
  fetch <- meta2$Length_km[i] * 1000 #sqrt((meta$SA.ha[i] * 1e5)/pi)*2 
  blen <- meta2$Length_km[i] * 1000 #sqrt((meta$SA.ha[i] * 1e5)/pi)*2
  bwid <- meta2$Width_km[i] * 1000 #sqrt((meta$SA.ha[i] * 1e5)/pi)*2
  
  input_yaml_multiple(file = ler_yaml, fetch, key1 = "model_parameters", key2 = "FLake", key3 = "fetch_lk")
  input_yaml_multiple(file = ler_yaml, blen, key1 = "model_parameters", key2 = "GLM", key3 = "bsn_len")
  input_yaml_multiple(file = ler_yaml, bwid, key1 = "model_parameters", key2 = "GLM", key3 = "bsn_wid")
  
  
  # set config
  config_file <- ler_yaml
  model <- c("GLM", "Simstrat")
  
  export_config(config_file = config_file, model = model)
  
  # run model ensemble
  run_ensemble(config_file = config_file, model = model, parallel = TRUE)
  
  ## not run - calibration 
  # cali_res <- cali_ensemble(config_file = config_file, num = 500, cmethod = "LHC",
  #                           parallel = TRUE, model = model, 
  #                           out_f = paste("cali", meta$lake.name[i], sep = "_"))
  # 
  # saveRDS(cali_res, paste("cali_", meta$lake.name[i], ".rds", sep = ""))
}

t2 <- Sys.time()
t2-t1
