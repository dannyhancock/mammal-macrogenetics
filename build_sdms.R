library(dplyr)
library(CoordinateCleaner)
library(dismo)
library(maptools)
library(raster)
library(biomod2)
library(remotes)
library(rgbif)
library(stringr)
source("functions.R")

### ENVIRONMENTAL DATA ###
path <- "maps/bioclim/bioclim_5k"
bioclim_files <- list.files(path, pattern='tif$',full.names=TRUE)
# create a rasterStack of bioclim predictor variables
bioclim_predictors <-  stack(bioclim_files)

hfp_path <- "maps/human_footprint/processed"
hfp_files <- list.files(hfp_path, pattern='tif$', full.names=TRUE)
hfp_predictors <- stack(hfp_files)

cop_path <- "maps/copernicus/processed"
cop_files <- list.files(cop_path, pattern='tif$', full.names=TRUE)
cop_predictors <- stack(cop_files)

### SPECIES OCCURRENCE DATA ###

# use the common name, see metadata.csv for mapping to scientific name
SPECIES <- 'american_pika'

gbif_download = readr::read_tsv(normalizePath(paste("./data/",SPECIES,"/",SPECIES,"_occurrence.csv", sep=""),
                              winslash = "/", mustWork = TRUE))

## CLEANING OCCURRENCE DATA ##
gbif_download <- gbif_download %>%
  setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
  filter(occurrencestatus  == "PRESENT") %>% # only keep presence data
  filter(!is.na(decimallongitude)) %>% # remove null coordinate occurrences, if any
  filter(!is.na(decimallatitude)) %>%
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
  filter(!establishmentmeans %in% c("MANAGED", "INTRODUCED", "INVASIVE", "NATURALISED")) %>%
  filter(year >= 1900) %>% # remove old records
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% # remove values with high uncertainty, missing values kept
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>% # remove values with high uncertainty, missing values kept
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% # known uncertainty values
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>% # remove points along prime meridian or equator
  distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>% # remove duplicates
  cc_cen(buffer = 2000, lon="decimallongitude", lat="decimallatitude") %>% # remove country centroids within 2km
  cc_cap(buffer = 2000, lon="decimallongitude", lat="decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon="decimallongitude", lat="decimallatitude") %>% # remove zoo and herbaria within 2km
  cc_sea(lon="decimallongitude", lat="decimallatitude") %>% # remove from ocean, as not marine species
  dplyr::rename(lon = decimallongitude, lat = decimallatitude) %>%
  glimpse()


## Add coordinates from study if few occurrences from GBIF ##
# sciureus_coords <- readr::read_csv(paste("data/",SPECIES,"/",SPECIES,"_indiv_locations.csv", sep=""))
# sciureus_coords$species <- 'Holochilus_sciureus'
# gbif_download <- rbind(unique(gbif_download[,c('lon', 'lat', 'species')]),sciureus_coords)
#
# dall_coords = readr::read_csv(paste("data/",SPECIES,"/",SPECIES,"_coords.csv", sep=""))
# dall_coords <- dall_coords[,c("longitude","latitude")]
# dall_coords <- cbind(dall_coords, "Ovis Dalli")
# colnames(dall_coords) <- c('lon','lat','species')
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],dall_coords)
# 
# prairie_coords <- readr::read_csv(paste("data/",SPECIES,"/","prairie_dogs_locations.csv", sep=""))
# prairie_coords <- prairie_coords[,c("longitude","latitude")]
# colnames(prairie_coords) <- c("lon","lat")
# prairie_coords$species <- 'Cynomys parvidens'
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],prairie_coords)
# 
# pygmy_coords <- readr::read_csv(paste("data/",SPECIES,"/","PopData_400.csv", sep=""))
# pygmy_coords <- unique(pygmy_coords[,c("longitude","latitude")])
# colnames(pygmy_coords) <- c("lon","lat")
# pygmy_coords$species <- 'Brachylagus idahoensis'
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],pygmy_coords)
# 
# snow_sheep_coords <- readr::read_csv(paste("data/",SPECIES,"/","Sampling_sites_coordinates.csv", sep=""))
# snow_sheep_coords <- unique(snow_sheep_coords[,c("Longitude","Latitude")])
# colnames(snow_sheep_coords) <- c("lon","lat")
# snow_sheep_coords$species <- 'Ovis nivicola'
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],snow_sheep_coords)
# 
## polar bear - remove wrapped coordinates ##
# gbif_download <- gbif_download %>%
#   filter(!lon > 0)
# 
# tas_coords <- readr::read_csv(paste("data/",SPECIES,"/","population_locations.csv", sep=""))
# tas_coords <- unique(tas_coords[,c("longitude","latitude")])
# colnames(tas_coords) <- c("lon","lat")
# tas_coords$species <- 'Sarcophilus harrisii'
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],tas_coords)
# 
# pebble_mouse_coords <- readr::read_csv(paste("data/",SPECIES,"/","SampleMetaData.csv", sep=""))
# pebble_mouse_coords <- pebble_mouse_coords[pebble_mouse_coords$species == "Pseudomys chapmani",]
# pebble_mouse_coords <- unique(pebble_mouse_coords[,c("lon","lat")])
# pebble_mouse_coords$species <- 'Pseudomys chapmani'
# gbif_download <- rbind(gbif_download[,c('lon', 'lat', 'species')],pebble_mouse_coords)

# quick visual check
data(wrld_simpl)
xlims <- c(min(gbif_download$lon)-10, max(gbif_download$lon)+10)
ylims <- c(min(gbif_download$lat)-10, max(gbif_download$lat)+10)
if (SPECIES == 'polar_bear') {
  xlims <- c(-200, -50)
  gbif_download <- filter(gbif_download, lon < -50)
}
plot(wrld_simpl, xlim=xlims, ylim=ylims, axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(gbif_download$lon, gbif_download$lat, col='orange', pch=20, cex=0.75)


### THINNING & OUTER TRAIN/TEST SPLIT ###
nrow(gbif_download)
acg <- gbif_download
coordinates(acg) <- ~lon+lat
crs(acg) <- crs(bioclim_predictors)
# create a raster layer with the extent of acg
r <- raster(acg)
# set the resolution of the cells to 0.5 degrees
res(r) <- 0.5
# expand the extent of the rasterlayer by 10 degrees
r <- extend(r, extent(r)+10)

# random sampling for training and test sets - 70/30 split
samp <- sample(nrow(acg), round(0.7 * nrow(acg)))
traindata <- acg[samp,]
traindata <- traindata@coords

testdata <- acg[-samp,]
testdata <- testdata@coords

# grid sampling to reduce sampling bias
acsel <- gridSample(traindata, r, n=2)
actest <- gridSample(testdata, r, n=1)

# plot the result
p <- rasterToPolygons(r)
plot(p, border='gray')
points(acg)
points(acsel, cex=1, col='red',pch='x')
points(actest, cex=1, col='blue',pch='x')

### MODELLING ###
## Prepare environment variables
# crop all predictors to extent of species + a bit
bioclim_cropped <- cropRasterStack(bioclim_predictors, r)
hfp_cropped <- resample(hfp_predictors, bioclim_cropped)
cop_cropped <- resample(cop_predictors, bioclim_cropped)

# landscape only
landscape_cropped <- addLayer(hfp_cropped, cop_cropped)

# all predictors
all_cropped <- addLayer(bioclim_cropped, hfp_cropped, cop_cropped)

# Build sdm per variable set
for (version in c('bioclim', 'landscape', 'all')){
  if (version == 'bioclim'){
    predictors_cropped <- bioclim_cropped
  } else if (version == 'landscape'){
    predictors_cropped <- landscape_cropped
  } else if (version == 'all'){
    predictors_cropped <- all_cropped
  }
  
  # produce presence vector for training and test sets
  myResp <- rep(1, nrow(acsel))
  myEval <- rep(1, nrow(actest))
  
  # first use BIOMOD to generate pseudoabsences for training data
  myBiomodPA_train <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = predictors_cropped,
                                           resp.xy = acsel,
                                           PA.nb.rep = 1,
                                           PA.nb.absences = length(myResp),
                                           PA.strategy = 'sre',
                                           PA.sre.quant = 0.1,
                                           resp.name = SPECIES)
  
  # grab pseudoabsence data from formatted dataset
  train_sp <- SpatialPointsDataFrame(coords = myBiomodPA_train@coord, 
                                     data = data.frame(species = myBiomodPA_train@data.species), 
                                     proj = CRS(proj4string(predictors_cropped)))
  
  # keep NAs for training data.
  # replace NAs with 0 for pseudobsences
  #train_sp$species[is.na(train_sp$species)] <- 0
  
  # first use BIOMOD to generate pseudoabsences for test data
  myBiomodPA_test <- BIOMOD_FormatingData(resp.var = myEval,
                                          expl.var = predictors_cropped,
                                          resp.xy = actest,
                                          resp.name = SPECIES,
                                          PA.nb.rep = 1,
                                          PA.nb.absences = length(myEval),
                                          PA.strategy = 'sre',
                                          PA.sre.quant = 0.1)
  # grab pseudoabsence data from formatted dataset
  test_sp <- SpatialPointsDataFrame(coords = myBiomodPA_test@coord, 
                                    data = data.frame(species = myBiomodPA_test@data.species), 
                                    proj = CRS(proj4string(predictors_cropped)))
  # replace NAs with 0 for pseudobsences
  test_sp$species[is.na(test_sp$species)] <- 0
  
  # check training datapoints
  plot(predictors_cropped[[1]])
  points(train_sp[!is.na(train_sp$species),], pch=19)
  points(train_sp[is.na(train_sp$species),], pch=24, col='red')
  
  # check test datapoints
  plot(predictors_cropped[[1]])
  points(test_sp[test_sp$species == 1,], pch=19)
  points(test_sp[test_sp$species == 0,], pch=24, col='red')
  
  # create background points equal to number of presence points
  background_points <- length(myResp)
  
  # now use BIOMOD to generate pseudoabsences for training data
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = predictors_cropped,
                                           resp.xy = acsel,
                                           PA.nb.rep = 5,
                                           PA.nb.absences = background_points,
                                           PA.strategy = 'sre',
                                           PA.sre.quant = 0.1,
                                           resp.name = SPECIES,
                                           eval.resp.var = test_sp$species,
                                           eval.expl.var = predictors_cropped,
                                           eval.resp.xy = test_sp@coords)
  myBiomodData

  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodData,
    models = c('GBM','RF'),
    OPT.strategy = "bigboss",
    #bm.options = myBiomodOption,
    CV.nb.rep = 5,
    CV.perc = 0.7, # split 70% of data for training in cross-validation
    prevalence = 0.5,
    var.import = 3,
    metric.eval = c('ROC','TSS','ACCURACY'),
    scale.models = TRUE,
    CV.do.full.models = FALSE,
    modeling.id = paste(SPECIES,"models",sep="_")
  )
  
  # get all models evaluation
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  myBiomodModelEval
  
  # plot eval boxplot
  bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
  
  # graph evalation scores by model
  bm_PlotEvalMean(
    myBiomodModelOut,
    dataset='evaluation',
    group.by = 'algo',
    metric.eval=c("ROC","TSS"),
    xlim=c(0,1),
    ylim=c(0,1)
  )
  
  # retrieve minimum score of best N models
  roc_df <- subset(myBiomodModelEval[c('full.name', 'algo', 'metric.eval', "calibration", 'validation', 'evaluation')], metric.eval=='ROC')
  
  # find top 5 models
  n_models <- 5
  roc_ordered <- roc_df[order(roc_df$evaluation, decreasing=TRUE),]
  top_5 <- roc_ordered[1:n_models,1]
  best_n <- min(tail(sort(roc_df$evaluation), n_models))-0.001
  
  # retrieve variable importance
  MyBiomodModelVarImp <- get_variables_importance(myBiomodModelOut)
  
  # get model names
  get_built_models(myBiomodModelOut)
  
  # variable importance for best individual model
  var_imp <- as.data.frame(MyBiomodModelVarImp)
  
  # Variable importance of best 5 models
  mean_variable_importance <- subset(var_imp, full.name %in% top_5) %>%
    group_by(expl.var) %>%
    summarise_at(vars(var.imp), list(name=mean))
  mean_variable_importance
  
  # Ensemble Modeling, choose models with highest ROC or best avg algorithm
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = top_5,
    em.by='all',
    metric.eval = c('ROC'),
    em.algo='EMwmean',
    committee.averaging = F
    )
  
  myBiomodEMEval <- get_evaluations(myBiomodEM)
  myBiomodEMEval
  
  # Ensemble by weighted mean ROC
  bestEns <- 'EMwmeanByROC'
  
  # Projections for average of best algorithm
  myBiomodAvgProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut,
    new.env = predictors_cropped,
    proj.name = 'projections',
    models.chosen = top_5,
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  
  # get projected map
  myCurrentProj <- get_predictions(myBiomodAvgProj)
  myCurrentProj
  
  # ensemble forecasting
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    bm.proj = myBiomodAvgProj,
    models.chosen = grep(bestEns,get_built_models(
      myBiomodEM), value=TRUE)
  )
  myBiomodEF
  # plot ensemble
  #plot(myBiomodEF)
  
  # conversion to geotiff for circuitscape
  biomodSPECIES <- str_replace_all(SPECIES, '_', '.')
  mygrd <- raster(paste(biomodSPECIES,"/proj_projections/proj_projections_", biomodSPECIES, "_ensemble.tif", sep=""))
  
  # Replace 0 suitability with minimum suitability (greater than 0)
  mygrd[is.na(mygrd)] <- min(mygrd[mygrd>0]) 
  mygrd[mygrd==0] <- min(mygrd[mygrd>0])
  
  plot(mygrd)
  mygrd
  myasc <- writeRaster(mygrd, paste("data/",SPECIES,"/circuitscape/",SPECIES,"_ensemble_",version,".asc", sep=""),format="ascii",overwrite=TRUE)
  
  # output results: models included in ensemble and AUC of ensemble.
  top_5_row <- roc_ordered[1:n_models,]
  top_5_algo <- top_5_row[,2]
  avg_roc_train <- mean(top_5_row$calibration)
  avg_roc_val <- mean(top_5_row$validation)
  avg_roc_eval <- mean(top_5_row$evaluation)
  
  ens <- subset(myBiomodEMEval, filtered.by == 'ROC' & algo == 'EMwmean')
  ens_roc_train <- ens$calibration
  ens_roc_val <- ens$validation
  ens_roc_eval <- ens$evaluation
  
  df_to_write = data.frame(SPECIES, version, 
                           top_5_algo[1],  top_5_algo[2], top_5_algo[3],
                           top_5_algo[4],  top_5_algo[5],
                           avg_roc_train, avg_roc_val, avg_roc_eval,
                           ens_roc_train, ens_roc_val, ens_roc_eval)
  
  write.table(df_to_write, file='sdm_results.csv', sep=',', append=TRUE,
              col.names=FALSE, row.names=FALSE)
  
  # Write variable importance to output
  imp_to_write <- data.frame(t(mean_variable_importance$name))
  colnames(imp_to_write) <- mean_variable_importance$expl.var
  meta <- data.frame(SPECIES, version)
  
  imp_to_write <- cbind(meta, imp_to_write)
  
  write.table(imp_to_write, file='variable_importance.csv', sep=',', append=TRUE,
              col.names=FALSE, row.names=FALSE)
  
  raster_check <- raster(paste("data/",SPECIES,"/circuitscape/",SPECIES,"_ensemble_",version,".asc", sep=""),)
  png(paste('sdm_images/',SPECIES, '_', version, '.png'))
  raster::plot(raster_check, main=paste(str_replace_all(SPECIES,'_',' '), '-', version, "variables"), cex.main=2, cex.axis=2, 
               axis.args=list(cex.axis=1.5))
  lines(wrld_simpl)
  dev.off()
}

