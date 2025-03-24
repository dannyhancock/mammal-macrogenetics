library(hierfstat)
library(vegan)
library(ecodist)
library(terra)
library(pegas)
library(stringr)
library(geosphere)
library(ResistanceGA)
source("functions.R")

# Load env data for IBE
bc_path <- "maps/bioclim/bioclim_5k"
bioclim_files <- list.files(bc_path, pattern='tif$',full.names=TRUE)
bioclim_unused <- list.files(paste(bc_path, "unused", sep="/"), full.names = TRUE)
all_bioclim_layers <- append(bioclim_files, bioclim_unused)
# create a rasterStack of predictor variables and crop to same extent as species occurences
bc_predictors <-  rast(all_bioclim_layers)
bc_predictors <- subset(bc_predictors, "wc2.1_2.5m_elev", negate=TRUE)
names(bc_predictors)

all_distances <- data.frame(
  Species = character(),
  genet_dist = numeric(),
  geog_dist = numeric(),
  env_dist = numeric(),
  biores_dist = numeric(),
  landres_dist = numeric(),
  allres_dist = numeric()
)

SCALE <- TRUE
# MOD <- "MRM"
MOD <- "LME"
EXTRA_ARG <- "_multi_optimized_scaled_AICc"

# Use same scaling as in MLPE.lmm function from ResistanceGA
scale_dist <- function(d) {
  m <- as.vector(d)
  m_scaled <- scale(m, center = TRUE, scale = TRUE)
  return(as.vector(m_scaled))
}

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
  geog_dist
  length(geog_dist)
  
  write.csv(as.data.frame(as.matrix(geog_dist)), paste0(path, SPECIES, "_geog_dist.csv"))

  ### ENVIRONMENTAL DISTANCE ##
  # extract environmental values at coordinates
  env_extract <- extract(bc_predictors, coords, method='simple')
  env_extract <- env_extract[, colnames(env_extract) %in% names(bc_predictors)]

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
  # If running on UNOPTIMIZED landscape based habitat suitability maps
  # landres <- read.table(paste(path, "circuitscape/", SPECIES, "_landscape_circuitscape_resistances.out", sep=""),
  #                         header=TRUE, row.names=1)
  # landres_dist <- as.dist(landres)
  # length(landres_dist)
  # landres_dist
  
  # If running on optimized landscape habitat suitability maps use this instead #
  landres <- read.table(paste0(path, "circuitscape/optimization/julia_parallel/Results/", SPECIES, "_ensemble_landscape_jlResistMat.csv"),
                        header=FALSE, row.names=NULL, sep = ",")
  landres_dist <- as.dist(landres)
  length(landres_dist)
  landres_dist
  
  ## all variables ##
  allres <- read.table(paste(path, "circuitscape/", SPECIES, "_all_circuitscape_resistances.out", sep=""),
                        header=TRUE, row.names=1)
  allres_dist <- as.dist(allres)
  length(allres_dist)
  allres_dist
  
  # Scale
  if (SCALE == TRUE) {
    geog_dist <- scale_dist(geog_dist)
    env_dist <- scale_dist(env_dist)
    biores_dist <- scale_dist(biores_dist)
    landres_dist <- scale_dist(landres_dist)
    allres_dist <- scale_dist(allres_dist)
  }
  
  # Create all possible combinations of models, e.g. genet_dist ~ geog_dist + env_dist
  distance_matrices <- list(
    geog_dist = geog_dist,
    env_dist = env_dist,
    biores_dist = biores_dist,
    landres_dist = landres_dist,
    allres_dist = allres_dist
  )
  
  # Generate all possible combinations of distance matrices
  all_combinations <- lapply(1:length(distance_matrices), function(i) {
    combn(names(distance_matrices), i, simplify = FALSE)
  }) |> unlist(recursive = FALSE)
  
  # Filter out combinations with more than one resistance distance
  valid_combinations <- Filter(function(combo) {
    sum(c("biores_dist", "landres_dist", "allres_dist") %in% combo) <= 1
  }, all_combinations)

  if (MOD == "LME") {
    ################################
    # Linear Mixed Effects Models ##
    ################################
    results_aic <- data.frame(SPECIES)
    
    for (distance in names(distance_matrices)) {
      
      # Model LME
      lme_model <- MLPE.lmm(distance_matrices[[distance]], lower_genet, REML = FALSE, scale=FALSE) # already scaled above
      
      aic_value <- AIC(logLik(lme_model))
      
      # Calculate AICc
      n <- nobs(lme_model)
      k <- attr(logLik(lme_model), "df")
      aic_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
      
      results_aic[[distance]] <- aic_value
    }
    
    # write results
    write.table(
      results_aic,
      paste0("results_LME/aic_results", EXTRA_ARG, ".csv"),
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(paste0("results_LME/aic_results", EXTRA_ARG, ".csv")),
      append = TRUE
    )
  }
  
  
  if (MOD == "MRM"){
    #########
    ## MRM ##
    #########
    results_r2 <- data.frame(SPECIES)
    results_r2_adj <- data.frame(SPECIES)
    results_p <- data.frame(SPECIES)
    results_F <- data.frame(SPECIES)
    
    for (combo in valid_combinations) {
      model_name <- paste(combo, collapse = " + ")
      
      formula <- as.formula(paste("genet_dist ~", paste(combo, collapse = " + ")))
      
      # Model MRM
      mod <- MRM(formula, nperm = 999, mrank = FALSE)
      
      # Calculate adjusted R^2
      n <- length(genet_dist)  # number of observations
      p <- length(combo)       # number of predictors
      r2 <- mod$r.squared[1]
      adj_r2 <- 1 - ((1 - r2) * (n - 1) / (n - p - 1))
      
      results_r2[[model_name]] <- mod$r.squared[1]
      results_r2_adj[[model_name]] <- adj_r2
      results_p[[model_name]] <- mod$r.squared[2]
      results_F[[model_name]] <- mod$F.test[1]
      
    }
    # write results
    write.table(results_r2, paste0('results_MRM/mrm_r2_results', EXTRA_ARG, '.csv'), sep = ",", row.names = FALSE, col.names = !file.exists(paste0('results_MRM/mrm_r2_results', EXTRA_ARG, '.csv')), append = TRUE)
    write.table(results_r2_adj, paste0('results_MRM/mrm_r2_adj_results', EXTRA_ARG, '.csv'), sep = ",", row.names = FALSE, col.names = !file.exists(paste0('results_MRM/mrm_r2_adj_results', EXTRA_ARG, '.csv')), append = TRUE)
    write.table(results_p, paste0('results_MRM/mrm_p_results', EXTRA_ARG, '.csv'), sep = ",", row.names = FALSE, col.names = !file.exists(paste0('results_MRM/mrm_p_results', EXTRA_ARG, '.csv')), append = TRUE)
    write.table(results_F, paste0('results_MRM/mrm_F_results', EXTRA_ARG, '.csv'), sep = ",", row.names = FALSE, col.names = !file.exists(paste0('results_MRM/mrm_F_results', EXTRA_ARG, '.csv')), append = TRUE)
    
  }

}