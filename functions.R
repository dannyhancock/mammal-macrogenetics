# Function to create a file path for a particular species for multiple uses
createpath <- function(species){
  path <- paste("data/", species, "/", sep="")
  return(path)
}

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

# Function for calculating ranks for models where delta AICc is within some limit
# are equivalent to the best model
adjusted_rank <- function(x, delta_limit=2) {
  # x is a vector of AIC values for a species
  best <- min(x)
  delta <- x - best
  # Models within delta_limit AIC units are tied for best:
  r <- rep(NA, length(x))
  best_idx <- which(delta <= delta_limit)
  r[best_idx] <- 1
  # For remaining models, order them by AIC and assign increasing rank starting at 2
  remaining_idx <- setdiff(seq_along(x), best_idx)
  if(length(remaining_idx) > 0) {
    ord <- order(x[remaining_idx])
    r[remaining_idx[ord]] <- seq(from = 2, length.out = length(remaining_idx))
  }
  return(r)
}

# Function to map models to labels for plotting
map_model_label <- function(model_str) {
  parts <- str_split(model_str, "_")[[1]]
  
  # Define mapping
  mapping <- list(
    geog = "IBD",
    env = "IBE",
    biores = "IBR[BIO]",
    landres = "IBR[LAND]",
    allres = "IBR[ALL]"
  )
  
  mapped <- sapply(parts, function(p) {
    if (!is.null(mapping[[p]])) mapping[[p]] else p
  })
  
  combined <- paste(mapped, collapse = " + ")
  
  parse(text = combined)
}

# Function for calculating confidence intervals for logistic regression
get_pred_ci <- function(model, newdata, level = 0.95) {
  # obtain predictions on the link scale along with standard errors
  pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  
  # calculate the critical value for the given confidence level (default 95%)
  crit_val <- qnorm((1 + level) / 2)
  
  # compute lower and upper bounds on the link scale
  fit_link <- pred$fit
  se_fit <- pred$se.fit
  lower_link <- fit_link - crit_val * se_fit
  upper_link <- fit_link + crit_val * se_fit
  
  # convert the estimates and bounds from the link (logit) scale to the response scale (probability)
  newdata$fit_prob <- plogis(fit_link)
  newdata$lower_prob <- plogis(lower_link)
  newdata$upper_prob <- plogis(upper_link)
  
  return(newdata)
}