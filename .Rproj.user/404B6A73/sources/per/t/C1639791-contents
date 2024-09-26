library(hierfstat)
library(vegan)
library(ecodist)
library(terra)
library(pegas)
library(stringr)
library(geosphere)

### IMPORTANT ### 
# ResistanceGA requires raster and sp however, sp interferes with ecodist which
# is required for MRM. The solution is to not load ResistanceGA and run MRM part first
# before MLPE.lmm part, then load ResistanceGA and run LME.

# library(ResistanceGA)

createpath <- function(species){
  path <- paste("data/", species, "/", sep="")
  return(path)
}

# Use common name from metadata.csv
SPECIES <- 'american_pika'

# Load env data for IBE
bc_path <- "maps/bioclim/bioclim_5k"
bioclim_files <- list.files(bc_path, pattern='tif$',full.names=TRUE)
bioclim_unused <- list.files(paste(bc_path, "unused", sep="/"), full.names = TRUE)
all_bioclim_layers <- append(bioclim_files, bioclim_unused)
# create a rasterStack of predictor variables and crop to same extent as species occurences
bc_predictors <-  rast(all_bioclim_layers)
bc_predictors <- subset(bc_predictors, "wc2.1_2.5m_elev", negate=TRUE)
names(bc_predictors)

dirs <- list.dirs("data", recursive = FALSE)
for (i in dirs){
  
  SPECIES <- str_split(i, "/", simplify=FALSE)[[1]][2]
  path <- createpath(SPECIES)
  
  ### GENETIC DISTANCE ###
  fst_csv <- read.csv(paste(path, SPECIES, '_fst_WC.csv', sep = ""), row.names=1)
  # replace negative values with 0
  fst_csv[fst_csv < 0] <- 0
  genet_dist <- as.dist(fst_csv, diag=FALSE)
  genet_dist
  length(genet_dist)

  lower_genet <- lower(as.matrix(genet_dist))

  ### GEOGRAPHIC DISTANCE ###
  coords <- read.table(paste0(path, "circuitscape/",SPECIES,"_focal_nodes.txt"), row.names = 1)
  geog_dist <- as.dist(geod(coords))
  geog_dist <- dist(coords, method="euclidean", diag=FALSE)
  geog_dist
  length(geog_dist)
  
  write.csv(as.data.frame(as.matrix(geog_dist)), paste0(path, SPECIES, "_geog_dist.csv"))

  ### ENVIRONMENTAL DISTANCE ##
  # extract environmental values at coordinates
  env_extract <- extract(bc_predictors, coords, method='simple')

  # PCA for feature extraction to determine environmental distance
  pca_mod <- rda(env_extract, scale=TRUE)
  summary(pca_mod)
  biplot(pca_mod)
  pca_scores <- scores(pca_mod,choices=c(1:2))$sites[,c(1:2)]
  pca_scores
  env_dist <- dist(pca_scores)
  env_dist
  length(env_dist)
  
  write.csv(as.data.frame(as.matrix(env_dist)), paste0(path, SPECIES, "_env_dist.csv"))

  pc2var <- summary(pca_mod)$cont$importance[3,"PC2"]
  
  ### RESISTANCE DISTANCE ###
  ## bioclim ##
  biores <- read.table(paste(path, "circuitscape/", SPECIES, "_bioclim_circuitscape_resistances.out", sep=""), 
                      header=TRUE, row.names=1)
  biores_dist <- as.dist(biores)
  length(biores_dist)
  biores_dist
  
  ## landscape ##
  landres <- read.table(paste(path, "circuitscape/", SPECIES, "_landscape_circuitscape_resistances.out", sep=""),
                          header=TRUE, row.names=1)
  landres_dist <- as.dist(landres)
  length(landres_dist)
  landres_dist
  
  ## all variables ##
  allres <- read.table(paste(path, "circuitscape/", SPECIES, "_all_circuitscape_resistances.out", sep=""),
                        header=TRUE, row.names=1)
  allres_dist <- as.dist(allres)
  length(allres_dist)
  allres_dist

  ## Linear Mixed Effects Models ##
  ## RUN MRM (below) FIRST ##

  # allres_lme <- MLPE.lmm(allres_dist,
  #                        lower_genet, REML=FALSE)
  # 
  # 
  # landres_lme <- MLPE.lmm(landres_dist,
  #                         lower_genet, REML=FALSE)
  # 
  # biores_lme <- MLPE.lmm(biores_dist,
  #                         lower_genet, REML=FALSE)
  # 
  # envres_lme <- MLPE.lmm(env_dist,
  #                         lower_genet, REML=FALSE)
  # 
  # geogres_lme <- MLPE.lmm(geog_dist,
  #                         lower_genet, REML=FALSE)
  # 
  # 
  # allres_lme
  # all_aic <- AIC(logLik(allres_lme))
  # landres_lme
  # land_aic <- AIC(logLik(landres_lme))
  # biores_lme
  # bio_aic <- AIC(logLik(biores_lme))
  # envres_lme
  # env_aic <- AIC(logLik(envres_lme))
  # geogres_lme
  # geog_aic <- AIC(logLik(geogres_lme))
  # 
  # df_to_write <- data.frame(SPECIES, all_aic, land_aic, bio_aic, env_aic, geog_aic, max_dist, pc2var)
  # 
  # write.table(df_to_write, "aic_results.csv", append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)

  # Multiple regression on distance matrices
  mod1 <- MRM(genet_dist ~ geog_dist, nperm=999, mrank=FALSE)
  mod1

  mod2 <- MRM(genet_dist ~ env_dist, nperm=999, mrank=FALSE)
  mod2

  mod3 <- MRM(genet_dist ~ biores_dist, nperm=999, mrank=FALSE)
  mod3

  mod4 <- MRM(genet_dist ~ landres_dist, nperm=999, mrank=FALSE)
  mod4

  mod5 <- MRM(genet_dist ~ allres_dist, nperm=999, mrank=FALSE)
  mod5

  mrm_r2_to_write <- data.frame(SPECIES, mod5$r.squared[1], mod4$r.squared[1], mod3$r.squared[1],
                                mod2$r.squared[1], mod1$r.squared[1])

  write.table(mrm_r2_to_write, 'mrm_r2_results.csv', append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)

  mrm_p_to_write <- data.frame(SPECIES, mod5$r.squared[2], mod4$r.squared[2], mod3$r.squared[2],
                               mod2$r.squared[2], mod1$r.squared[2])

  write.table(mrm_p_to_write, 'mrm_p_results.csv', append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)

  mrm_F_to_write <- data.frame(SPECIES, mod5$F.test[1], mod4$F.test[1], mod3$F.test[1],
                               mod2$F.test[1], mod1$F.test[1])

  write.table(mrm_F_to_write, 'mrm_F_results.csv', append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)

}