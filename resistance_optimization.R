library(raster)
library(ResistanceGA)
library(sp)
library(tidyverse)
library(maptools)

data(wrld_simpl)

JULIA_HOME="C:/Users/user/AppData/Local/Programs/Julia-1.10.2/bin"
JuliaCall::julia_setup(JULIA_HOME)

dirs <- list.dirs("data", recursive = FALSE)
times_GA <- list()  
for (species_dir in dirs) {
  SPECIES <- basename(species_dir)
  
  cat("Processing species: ", SPECIES, "\n")
  
  Results_dir <- paste0("data/",SPECIES, "/circuitscape/optimization/julia_parallel/")
  # Create the directory if it doesn't exist
  if (!dir.exists(Results_dir)) {
    dir.create(Results_dir, recursive = TRUE)
  }
  
  # Genetic distance matrix
  gen_file <- file.path(species_dir, paste0(SPECIES, "_fst_WC.csv"))
  genetic_distances <- as.dist(read.csv(gen_file, row.names = 1))
  Fst <- as.matrix(genetic_distances)
  
  # Spatial coordinates
  CS_Point.file <- paste0("data/", SPECIES, "/circuitscape/", SPECIES, "_focal_nodes.txt")
  coordinates <- read.table(CS_Point.file, header=F, row.names=1)
  sample.locales <- SpatialPoints(coordinates)
  
  # Prepare the Genetic Algorithm inputs.
  raster_file <- file.path(species_dir, "circuitscape", paste0(SPECIES, "_ensemble_landscape", ".asc"))
  r <- raster(raster_file)
  
  e <- extent(sample.locales)
  e <- e + 15
  r_crop <- crop(r, e)
  ncells <- ncell(r_crop)
  cat("Number of cells: ", ncells, "\n")
  
  plot(r_crop)
  points(sample.locales, cex=0.9, pch=20, col="blue")
  
  # Prepare inputs
  GA.inputs <- GA.prep(ASCII.dir = r_crop,
                       Results.dir = Results_dir,
                       maxiter = 500,
                       run = 25,
                       quiet = FALSE,
                       select.trans = list("A"),
                       parallel=6,
                       seed=1337)
  
  jl.inputs <- jl.prep(n.Pops = length(sample.locales),
                       response = lower(Fst),
                       CS_Point.File = sample.locales,
                       JULIA_HOME = JULIA_HOME)
  
  # Optimize
  start_time <- Sys.time()
  
  jl.optim <- SS_optim(jl.inputs = jl.inputs,
                       GA.inputs = GA.inputs)
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
  
  time_taken_mins <- as.numeric(time_taken, units = "mins")
  
  # Write timing of optimization to csv
  write.table(data.frame(SPECIES = SPECIES,
                         ncells = ncells,
                         time_taken_mins = time_taken_mins),
              "optimization_times_GA.csv",
              sep=",",
              append = TRUE,
              row.names=FALSE,
              col.names = !file.exists("optimization_times_GA.csv"))
  
  optim_resistance <- raster(paste0(Results_dir, "Results/", SPECIES, "_ensemble_landscape.asc"))
  plot(optim_resistance)
}
