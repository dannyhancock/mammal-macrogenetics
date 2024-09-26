library(raster)

# function to project and aggregate rasters to same as bioclim rasters
resizeRasterStack <- function(source_stack, target_stack, agg_factor=4) {
  print(names(source_stack))
  if (identical(crs(source_stack), crs(target_stack))) {
    print("Identical CRS, Aggregating")
    resolution <- res(target_stack)
    fact <- resolution/res(source_stack)
    agg_raster <- aggregate(source_stack, fact, fun='mean')
  } else {
    print("Raster Stack Requires Reprojecting")
    gc()
    agg_raster <- aggregate(source_stack, fact=agg_factor) # this line is to make the projection manageable with small amounts of RAM but may need changing
    agg_raster <- projectRaster(agg_raster, res=res(target_stack[[1]]), crs=crs(target_stack[[1]]))
  }
  return(agg_raster)
}

# function for cropping raster stack with output of rasterstack 
# (rather than rasterbrick) for input into biomod2
cropRasterStack <- function(stack, e) {
  cropped_predictors <- stack()
  for (i in 1:length(names(stack))) {
    cropped_layer <- crop(stack[[i]], extent(e))
    cropped_predictors <- addLayer(cropped_predictors, cropped_layer)
  }
  return(cropped_predictors)
}


##################
## Bioclim Maps ##
##################
# Other maps will be projected and aggregated to match Bioclim maps:
# WGS84 projection 5km resolution
bc_path <- "maps/bioclim/bioclim_5k"
bc_files <- list.files(bc_path, pattern='tif$', full.names=TRUE)
bc_predictors <- stack(bc_files)

##########################
## Human Footprint Maps ##
##########################
raw_hf_path <- "maps/human_footprint/raw"
raw_hf_files <- list.files(raw_hf_path, pattern='tif$', full.names=TRUE)
raw_hf_predictors <- stack(raw_hf_files)

# Reproject Human Footprint Maps to WGS84
hfp_predictors <- projectRaster(raw_hf_predictors, bc_predictors[[1]])

#####################
## Copernicus Maps ##
#####################
raw_cop_path <- "maps/copernicus/raw"
raw_cop_files <- list.files(raw_cop_path, pattern='tif$', full.names=TRUE)
raw_cop_predictors <- stack(raw_cop_files)
raw_cop_predictors

cop_predictors <- resizeRasterStack(raw_cop_predictors, bc_predictors)

##################
## Collinearity ##
##################
# Read in all species GBIF coordinates
dirs <- list.dirs(path = "data/", full.names = FALSE, recursive = FALSE)
allsp <- data.frame()
for (dir in dirs) {
  sp <- read.csv(paste("data", dir, paste0(dir, "_occurrence.csv"), sep="/"))
  allsp <- rbind(allsp, sp)
}

e <- extent(x=c(ceiling(min(allsp$lon)), ceiling(max(allsp$lon))), y=c(ceiling(min(allsp$lat)),ceiling(max(allsp$lat))))+5
# crop all predictor stacks to same extent to allow stacking
cropped_bioclim <- cropRasterStack(bc_predictors, e)
cropped_hfp <- resample(hfp_predictors, cropped_bioclim)
cropped_cop <- resample(cop_predictors, cropped_bioclim)

# stack predictors variables
cropped_predictors <- addLayer(cropped_bioclim, cropped_hfp, cropped_cop)

# extract coordinates and look for correlations
sampl <- cbind(allsp$lon, allsp$lat)
sampled_preds <- extract(cropped_predictors, sampl, 'bilinear')

cor_all <- cor(sampled_preds, method = 'pearson', use='complete.obs')
cor_all > 0.7
cor_all < -0.7

# keep most biologically interpretable uncorrelated layers
bioclim_predictors <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_2', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_8',
                        'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_18')
uncorrelated_bc_predictors <- subset(bc_predictors, bioclim_predictors)

uncorrelated_hfp_predictors <- subset(hfp_predictors, c('Built2009','croplands2005','Pasture2009','Railways', 'Roads'))

uncorrelated_cop_predictors <- subset(cop_predictors, c('grass', 'shrub', 'tree'))

# Move the correlated bioclim predictors into a separate folder
bioclim_dir <- 'maps/bioclim/bioclim_5k'
unused_dir = 'maps/bioclim/bioclim_5k/unused'
files_in_bioclim <- list.files(bioclim_dir, pattern='tif$', full.names=F)
# append .tif to the predictor names
bioclim_predictors_tif <- paste0(bioclim_predictors, ".tif")
files_to_move <- setdiff(files_in_bioclim, bioclim_predictors_tif)
for (file in files_to_move) {
  file_path <- file.path(bioclim_dir, file)
  dest_path <- file.path(unused_dir, file)
  
  if (!dir.exists(unused_dir)) {
    dir.create(unused_dir)
  }
  
  # Move the file
  file.rename(file_path, dest_path)
  
  cat(paste("Moved:", file, "to unused folder.\n"))
}

# Write the processed human footprint predictors to disk
output_path = './maps/human_footprint/processed/'
for (i in 1:nlayers(uncorrelated_hfp_predictors)) {
  layer_name <- paste(output_path, names(uncorrelated_hfp_predictors)[i], '.tif', sep = "")
  writeRaster(uncorrelated_hfp_predictors[[i]], filename = layer_name, format = "GTiff", overwrite = TRUE)
}

# and write the Copernicus predictors
output_path = './maps/copernicus/processed'
for (i in 1:nlayers(uncorrelated_cop_predictors)) {
  layer_name <- paste(output_path, names(uncorrelated_cop_predictors)[i], '.tif', sep = "")
  print(layer_name)
  raster::writeRaster(uncorrelated_cop_predictors[[i]], filename = layer_name, format = "GTiff", overwrite = TRUE)
}
