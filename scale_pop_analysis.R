# This sript is to conduct logistic regression to model the effect of
# the number of populations and the spatial scale on the likelihood of
# (IBD, IBE, IBR) being chosen as, or included in, the best model of genetic
# divergence in mammals

library(readxl)
library(tidyr)
library(nnet)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
source("functions.R")

metadata_sheet <- read.csv("metadata.csv")
metadata <- metadata_sheet[,c('species', 'populations', 'max_dist')]

##############
## R2 RANKS ##
##############
# extra_arg <- "_multi_unoptimized_scaled"
extra_arg <- "_multi_optimized_scaled"
r_scores <- read.csv(paste0('mrm_r2_adj_results', extra_arg, '.csv'))

# quick checks
identical(nrow(metadata), nrow(r_scores))
identical(metadata$species, r_scores$species)

cols_to_keep <- c("geog_dist", "env_dist", "landres_dist",
                  "geog_dist...env_dist", "geog_dist...landres_dist", "env_dist...landres_dist",
                  "geog_dist...env_dist...landres_dist")

mod_order <- c("geog", "env", "landres",
               "geog_env", "geog_landres", "env_landres",
               "geog_env_landres")

# remove allres and biores and apply ranking
r_ranks <- t(apply(-r_scores[,cols_to_keep], 1, rank))
colnames(r_ranks) <- mod_order
r_ranks_df <- data.frame(species = r_scores$SPECIES, r_ranks)

# merge with metadata to combine no. populations and maximum distance between pops
r_df <- merge(r_ranks_df, metadata, by = 'species')

# find the top ranked driver for each species
r_df$top_ranked <- apply(r_df[, mod_order], 1, function(x) mod_order[which.min(x)])
r_df$top_ranked <- factor(r_df$top_ranked, levels = mod_order)

r_df$has_geog    <- ifelse(grepl("geog", r_df$top_ranked), 1, 0)
r_df$has_env     <- ifelse(grepl("env", r_df$top_ranked), 1, 0)
r_df$has_landres <- ifelse(grepl("landres", r_df$top_ranked), 1, 0)

colSums(r_df[,c('has_geog', 'has_env', 'has_landres')])

### No. populations ###
# Drop outlier as logistic regression is sensitive to extreme values
pop_r_df <- r_df[r_df$species != 'bank_vole',]

r_pop_model_geog <- glm(has_geog ~ populations, data = pop_r_df, family = binomial)
r_pop_model_env  <- glm(has_env  ~ populations, data = pop_r_df, family = binomial)
r_pop_model_landres <- glm(has_landres ~ populations, data = pop_r_df, family = binomial)

summary(r_pop_model_geog)
summary(r_pop_model_env)
summary(r_pop_model_landres)

r_pop_seq <- seq(min(pop_r_df$populations), max(pop_r_df$populations), length.out = 120)

# Create a new data frame for predictions:
r_pop_new_data <- data.frame(populations = r_pop_seq)

# Predict probabilities for each predictor:
r_pop_pred_geog <- get_pred_ci(r_pop_model_geog, r_pop_new_data)
r_pop_pred_geog$Model <- "geog"
r_pop_pred_env  <- get_pred_ci(r_pop_model_env, r_pop_new_data)
r_pop_pred_env$Model <- "env"
r_pop_pred_landres <- get_pred_ci(r_pop_model_landres, r_pop_new_data)
r_pop_pred_landres$Model <- "landres"

r_pop_pred_all <- rbind(r_pop_pred_geog, r_pop_pred_env, r_pop_pred_landres)
r_pop_pred_all$Model <- factor(r_pop_pred_all$Model, levels = c("geog", "env", "landres"))

# Plot all models on the same plot
r_pop_plot <- ggplot(r_pop_pred_all, aes(x = populations, y = fit_prob, color = Model, fill = Model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA) +
  ylim(0, 1) +
  labs(x = "Number of Populations", y = "Probability", title = "MRM") +
  scale_color_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  scale_fill_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
r_pop_plot

### Spatial Scale ###
# Fit the logistic regression model
r_dist_model_geog <- glm(has_geog ~ max_dist, data = r_df, family = binomial)
r_dist_model_env  <- glm(has_env  ~ max_dist, data = r_df, family = binomial)
r_dist_model_landres <- glm(has_landres ~ max_dist, data = r_df, family = binomial)

summary(r_dist_model_geog)
summary(r_dist_model_env)
summary(r_dist_model_landres)

r_dist_seq <- seq(min(r_df$max_dist), max(r_df$max_dist), length.out = 120)

# Create a new data frame for predictions:
r_dist_new_data <- data.frame(max_dist = r_dist_seq)

# Predict probabilities for each predictor:
r_dist_pred_geog <- get_pred_ci(r_dist_model_geog, r_dist_new_data)
r_dist_pred_geog$Model <- "geog"
r_dist_pred_env  <- get_pred_ci(r_dist_model_env, r_dist_new_data)
r_dist_pred_env$Model <- "env"
r_dist_pred_landres <- get_pred_ci(r_dist_model_landres, r_dist_new_data)
r_dist_pred_landres$Model <- "landres"

r_dist_pred_all <- rbind(r_dist_pred_geog, r_dist_pred_env, r_dist_pred_landres)
r_dist_pred_all$Model <- factor(r_dist_pred_all$Model, levels = c("geog", "env", "landres"))

# Plot predicted probabilities vs. number of populations:
r_dist_plot <- ggplot(r_dist_pred_all, aes(x = max_dist, y = fit_prob, color = Model, fill=Model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA) +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle("MRM")+
  labs(x = "Maximum Distance (Km)",
       y = "Probability") +
  scale_color_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label))+
  scale_fill_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  theme(plot.title = element_text(hjust = 0.5))
r_dist_plot

########################
# Repeat for AIC ranks #
########################
extra_arg <- "_multi_optimized_scaled_AICc"
aic_scores <- read.csv(paste0('aic_results', extra_arg, '.csv'))

# quick checks
identical(nrow(metadata), nrow(aic_scores))
identical(metadata$species, aic_scores$species)

# remove allres and biores and apply ranking
cols_to_keep <- c("geog_dist", "env_dist", "landres_dist")
aic_ranks <- t(apply(aic_scores[,cols_to_keep], 1, function(x) adjusted_rank(x, 7)))
colnames(aic_ranks) <- cols_to_keep
aic_ranks_df <- data.frame(species = aic_scores$SPECIES, aic_ranks)

# merge with metadata to combine no. populations and maximum distance between pops
aic_df <- merge(aic_ranks_df, metadata, by = 'species')

aic_df$has_geog    <- ifelse(aic_df$geog_dist    == 1, 1, 0)
aic_df$has_env     <- ifelse(aic_df$env_dist     == 1, 1, 0)
aic_df$has_landres <- ifelse(aic_df$landres_dist == 1, 1, 0)

### No. populations ###
# Drop outlier as multinomial logistic regression is sensitive to extreme values
aic_pop_df <- aic_df[aic_df$species != 'bank_vole',]

aic_pop_model_geog    <- glm(has_geog ~ populations, data = aic_pop_df, family = binomial)
aic_pop_model_env     <- glm(has_env ~ populations, data = aic_pop_df, family = binomial)
aic_pop_model_landres <- glm(has_landres ~ populations, data = aic_pop_df, family = binomial)

# Inspect summaries
summary(aic_pop_model_geog)
summary(aic_pop_model_env)
summary(aic_pop_model_landres)

aic_pop_pop_seq <- seq(min(aic_pop_df$populations), max(aic_pop_df$populations), length.out = 120)
aic_pop_new_data <- data.frame(populations = aic_pop_pop_seq)

# Predict probabilities for each model
# Predict probabilities for each predictor:
aic_pop_pred_geog <- get_pred_ci(aic_pop_model_geog, aic_pop_new_data)
aic_pop_pred_geog$Model <- "geog"
aic_pop_pred_env  <- get_pred_ci(aic_pop_model_env, aic_pop_new_data)
aic_pop_pred_env$Model <- "env"
aic_pop_pred_landres <- get_pred_ci(aic_pop_model_landres, aic_pop_new_data)
aic_pop_pred_landres$Model <- "landres"

aic_pop_pred_all <- rbind(aic_pop_pred_geog, aic_pop_pred_env, aic_pop_pred_landres)
aic_pop_pred_all$Model <- factor(aic_pop_pred_all$Model, levels = c("geog", "env", "landres"))

# Plot the predicted probabilities vs. number of populations
aic_pop_plot <- ggplot(aic_pop_pred_all, aes(x = populations, y = fit_prob, color = Model, fill=Model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA) +
  ylim(0, 1) +
  labs(x = "Number of Populations",
       y = "Probability") +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~7)) +
  scale_color_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  scale_fill_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
aic_pop_plot

## Spatial scale ##
aic_dist_model_geog    <- glm(has_geog ~ max_dist, data = aic_df, family = binomial)
aic_dist_model_env     <- glm(has_env ~ max_dist, data = aic_df, family = binomial)
aic_dist_model_landres <- glm(has_landres ~ max_dist, data = aic_df, family = binomial)

# Inspect summaries
summary(aic_dist_model_geog)
summary(aic_dist_model_env)
summary(aic_dist_model_landres)

aic_dist_dist_seq <- seq(min(aic_df$max_dist), max(aic_df$max_dist), length.out = 120)
aic_dist_new_data <- data.frame(max_dist = aic_dist_dist_seq)

# Predict probabilities for each model
aic_dist_pred_geog <- get_pred_ci(aic_dist_model_geog, aic_dist_new_data)
aic_dist_pred_geog$Model <- "geog"
aic_dist_pred_env  <- get_pred_ci(aic_dist_model_env, aic_dist_new_data)
aic_dist_pred_env$Model <- "env"
aic_dist_pred_landres <- get_pred_ci(aic_dist_model_landres, aic_dist_new_data)
aic_dist_pred_landres$Model <- "landres"

aic_dist_pred_all <- rbind(aic_dist_pred_geog, aic_dist_pred_env, aic_dist_pred_landres)
aic_dist_pred_all$Model <- factor(aic_dist_pred_all$Model, levels = c("geog", "env", "landres"))

# Plot the predicted probabilities vs. the maximum distance
aic_dist_plot <- ggplot(aic_dist_pred_all, aes(x = max_dist, y = fit_prob, color = Model, fill=Model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA) +
  ylim(0, 1) +
  labs(x = "Maximum Distance (Km)",
       y = "Probability") +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~7)) +
  scale_color_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label))+
  scale_fill_brewer(palette = "Set2", labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())+
  annotate("text", x = 6200, y = 0.97, label = "*", size = 10, color = "black", fontface = "bold")
aic_dist_plot

# Arrange the plots
aic_dist_plot <- aic_dist_plot + theme(axis.title.y = element_blank(),
                                       axis.title.x = element_blank())
aic_pop_plot <- aic_pop_plot + theme(axis.title.x = element_blank())

combined_plot <- ggarrange(
  aic_pop_plot, aic_dist_plot, r_pop_plot, r_dist_plot, 
  ncol = 2, nrow = 2, 
  common.legend = TRUE, legend = "right",
  labels = c("A", "B", "C", "D"),
  label.x = 0.05
)

print(combined_plot)

ggsave("plots/dist_pop_plots.svg", plot = combined_plot, width = 10, height = 8)
