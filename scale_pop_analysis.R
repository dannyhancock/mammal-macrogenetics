# This sript is to conduct multinomial logistic regression to model the effect of
# the number of populations on the likelihood of a driver of genetic differentiation
# (IBD, IBE, IBR) being chosen as the best model.

library(readxl)
library(tidyr)
library(nnet)
library(tidyr)
library(ggplot2)
library(ggpubr)

calculate_p_values <- function(model){
  mod_summary <- summary(model)
  
  # Calculate z-values
  z_values <- mod_summary$coefficients / mod_summary$standard.errors
  
  # Calculate p-values with 2-tailed z-test
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  return(p_values)
}

plot_multinom <- function(pred_probs_long, x_var, sig, sig_size=5) {
  x_label <- if (x_var == "populations") {
    "Number of Populations"
  } else if (x_var == "max_dist") {
    "Euclidean Distance (degrees)"
  } else {
    x_var
  }
  
  # Plot
  p <- ggplot(pred_probs_long, aes_string(x = x_var, y = "Probability", color = "Category")) +
    geom_line(size = 1) +
    ylim(0, 1) +
    scale_color_brewer(palette = "Dark2", name = "Driver") +
    labs(
      x = x_label, 
      y = "Probability"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, face = "bold")
    ) +
    annotate("text", x = Inf, y = Inf, label = sig, hjust = 1.5, vjust = 2.5, size = sig_size, fontface = "bold")
  
  return (p)
}

metadata_sheet <- read_excel("metadata.csv")
metadata <- metadata_sheet[,c('species', 'populations', 'max_dist')]

##############
## R2 RANKS ##
##############
r_scores <- read.csv('mrm_r2_results.csv')

# quick checks
identical(nrow(metadata), nrow(r_scores))
identical(metadata$species, r_scores$species)

# remove allres and biores and apply ranking
r_ranks <- t(apply(-r_scores[,c('landres', 'env', 'geog')], 1, rank))
r_ranks_df <- data.frame(species = r_scores$species, r_ranks)

# merge with metadata to combine no. populations and maximum distance between pops
r_df <- merge(r_ranks_df, metadata, by = 'species')

# find the top ranked driver for each species
ranked_columns <- c("landres", "env", "geog")
r_df$top_ranked <- apply(r_df[, ranked_columns], 1, function(x) ranked_columns[which.min(x)])
r_df$top_ranked <- factor(r_df$top_ranked, levels = c("landres", "env", "geog"))


### No. populations ###
# Drop outlier as multinomial logistic regression is sensitive to extreme values
pop_r_df <- r_df[r_df$species != 'bank_vole',]

pop_mod <- nnet::multinom(top_ranked ~ populations, data = pop_r_df)
summary(pop_mod)
pop_pv <- calculate_p_values(pop_mod)
print("p values for effect of population (MRM ranks):")
pop_pv

# Calculating Odds Ratios
exp(coef(pop_mod))

# Plot the results
# Generate a sequence of values for populations
pop_seq <- seq(min(pop_r_df$populations), max(pop_r_df$populations), length.out = 120)

# Create a new data frame for predictions
new_data <- data.frame(populations = pop_seq)

# Predict probabilities
pred_probs <- predict(pop_mod, new_data, type = "probs")

# Convert predictions to long format for plotting

pred_probs_long <- gather(data.frame(populations = pop_seq, pred_probs), 
                          key = "Category", value = "Probability", -populations)

pred_probs_long$Category <- factor(pred_probs_long$Category,
                                        levels = c("geog", "env", "landres"),
                                        labels = c("IBD", "IBE", "IBR"))

# Plot
mrm_popplot <- plot_multinom(pred_probs_long, "populations", "n.s.")
mrm_popplot

### Spatial Scale ###
# Fit the multinomial logistic regression model
dist_mod <- nnet::multinom(top_ranked ~ as.numeric(max_dist), data = r_df)
summary(dist_mod)
dist_mod

spatial_pv <- calculate_p_values(dist_mod)
print("p values effect of spatial scale (MRM rankings):")
spatial_pv

# Calculating Odds Ratios
exp(coef(dist_mod))

pp <- fitted(dist_mod)
pp

# Plot the results
# Generate a sequence of values for populations
dist_pop_seq <- seq(min(r_df$max_dist), max(r_df$max_dist), length.out = 120)

# Create a new data frame for predictions
dist_new_data <- data.frame(max_dist = dist_pop_seq)

# Predict probabilities
dist_pred_probs <- predict(dist_mod, dist_new_data, type = "probs")

# Convert predictions to long format for plotting
dist_pred_probs_long <- gather(data.frame(max_dist = dist_pop_seq, dist_pred_probs), 
                          key = "Category", value = "Probability", -max_dist)

dist_pred_probs_long$Category <- factor(dist_pred_probs_long$Category,
                                        levels = c("geog", "env", "landres"),
                                        labels = c("IBD", "IBE", "IBR"))


mrm_distplot <- plot_multinom(dist_pred_probs_long, "max_dist", "*  ", sig_size=10)
mrm_distplot


########################
# Repeat for AIC ranks #
########################
aic_scores <- read.csv('aic_results.csv')

# quick checks
identical(nrow(metadata), nrow(aic_scores))
identical(metadata$species, aic_scores$species)

# remove allres and biores and apply ranking
aic_ranks <- t(apply(aic_scores[,c('landres', 'env', 'geog')], 1, rank))
aic_ranks_df <- data.frame(species = aic_scores$species, aic_ranks)

# merge with metadata to combine no. populations and maximum distance between pops
aic_df <- merge(aic_ranks_df, metadata, by = 'species')

aic_df <- aic_df[aic_df$species != 'sugar_glider',]
# find the top ranked driver for each species
aic_df$top_ranked <- apply(aic_df[, ranked_columns], 1, function(x) ranked_columns[which.min(x)])
aic_df$top_ranked <- factor(aic_df$top_ranked, levels = c("landres", "env", "geog"))


### No. populations ###
# Drop outlier as multinomial logistic regression is sensitive to extreme values
aic_df <- aic_df[aic_df$species != 'bank_vole',]

aic_pop_mod <- nnet::multinom(top_ranked ~ populations, data = aic_df)
summary(aic_pop_mod)
aic_pop_pv <- calculate_p_values(aic_pop_mod)
aic_pop_pv

# Calculating Odds Ratios
exp(coef(aic_pop_mod))

# Plot the results
# Generate a sequence of values for populations
aic_pop_seq <- seq(min(aic_df$populations), max(aic_df$populations), length.out = 120)

# Create a new data frame for predictions
aic_new_data <- data.frame(populations = aic_pop_seq)

# Predict probabilities
aic_pred_probs <- predict(aic_pop_mod, aic_new_data, type = "probs")

# Convert predictions to long format for plotting

aic_pop_pred_probs_long <- gather(data.frame(populations = aic_pop_seq, aic_pred_probs), 
                          key = "Category", value = "Probability", -populations)
aic_pop_pred_probs_long$Category <- factor(aic_pop_pred_probs_long$Category,
                                            levels = c("geog", "env", "landres"),
                                            labels = c("IBD", "IBE", "IBR"))

aic_popplot <- plot_multinom(aic_pop_pred_probs_long, "populations", "n.s.")
aic_popplot

### Spatial Scale ###
# Fit the multinomial logistic regression model
aic_dist_mod <- nnet::multinom(top_ranked ~ as.numeric(max_dist), data = aic_df)
summary(aic_dist_mod)

aic_spatial_pv <- calculate_p_values(aic_dist_mod)
aic_spatial_pv

# Calculating Odds Ratios
exp(coef(aic_dist_mod))

pp <- fitted(aic_dist_mod)
pp

# Plot the results
# Generate a sequence of values for populations
aic_dist_seq <- seq(min(aic_df$max_dist), max(aic_df$max_dist), length.out = 105)

# Create a new data frame for predictions
aic_dist_new_data <- data.frame(max_dist = aic_dist_seq)

# Predict probabilities
aic_dist_pred_probs <- predict(aic_dist_mod, aic_dist_new_data, type = "probs")

# Convert predictions to long format for plotting
aic_dist_pred_probs_long <- gather(data.frame(max_dist = aic_dist_seq, aic_dist_pred_probs), 
                               key = "Category", value = "Probability", -max_dist)

aic_dist_pred_probs_long$Category <- factor(aic_dist_pred_probs_long$Category,
                                        levels = c("geog", "env", "landres"),
                                        labels = c("IBD", "IBE", "IBR"))

aic_distplot <- plot_multinom(aic_dist_pred_probs_long, "max_dist", "n.s.")
aic_distplot

# Arrange the plots
mrm_distplot <- mrm_distplot + ggtitle("MRM")
aic_distplot <- aic_distplot + ggtitle("LME") + theme(axis.title.y = element_blank())
aic_popplot <- aic_popplot + theme(axis.title.y = element_blank(), plot.margin=margin(20,0,0,0))
mrm_popplot <- mrm_popplot + theme(plot.margin=margin(20,0,0,0))

combined_plot <- ggarrange(
  mrm_distplot, aic_distplot, mrm_popplot, aic_popplot, 
  ncol = 2, nrow = 2, 
  common.legend = TRUE, legend = "right",
  labels = c("A", "B", "C", "D")
)

print(combined_plot)

ggsave("plots/driver_multinom.svg", plot = combined_plot, width = 10, height = 8)
