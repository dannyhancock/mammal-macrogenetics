library(ggplot2)
library(readxl)
library(plyr)
library(ggpubr)
library(data.table)
library(gridExtra)
library(maptools)
library(cowplot)
source("functions.R")

#########################
### MODELLING RESULTS ###
#########################

## Resistance Only - which (unoptimzied) resistance model best explains genetic divergence ##
extra_arg <- "_multi_unoptimized_scaled_AICc"
modresults <- read.csv(paste0('aic_results', extra_arg, '.csv'))

aic_resonly <- modresults[,c('biores_dist','landres_dist','allres_dist')]
aic_mat <- as.matrix(aic_resonly[, c("biores_dist", "landres_dist", "allres_dist")])

# first for delta 2
aic_resonly_adjusted_ranks_d2 <- t(apply(aic_mat, 1, function(x) adjusted_rank(x, delta_threshold=2)))
aic_resonly_ranksums_adj_d2 <- colSums(aic_resonly_adjusted_ranks_d2 ==1)

first_ranks_d2 <- colSums(aic_resonly_adjusted_ranks_d2 == 1)

# Convert to a data frame for ggplot
df_first_d2 <- data.frame(
  Model = factor(c("biores", "landres", "allres"), levels = c("biores", "landres", "allres")),
  Count = as.numeric(first_ranks_d2)
)

# Plot the counts
aic_resonly_plot_d2 <- ggplot(df_first_d2, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Blues") +
  ylim(0, 40) +
  xlab("Model") +
  ylab("Count of Datasets with Lowest AICc") +
  theme_minimal() +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~2)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
aic_resonly_plot_d2

# repeat for delta 7
aic_resonly_adjusted_ranks_d7 <- t(apply(aic_mat, 1, function(x) adjusted_rank(x, delta_threshold=7)))
aic_resonly_ranksums_adj_d7 <- colSums(aic_resonly_adjusted_ranks_d7 ==1)

first_ranks_d7 <- colSums(aic_resonly_adjusted_ranks_d7 == 1)

# Convert to a data frame for ggplot
df_first_d7 <- data.frame(
  Model = factor(c("biores", "landres", "allres"), levels = c("biores", "landres", "allres")),
  Count = as.numeric(first_ranks_d7)
)

aic_resonly_plot_d7 <- ggplot(df_first_d7, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Blues", name = "IBR Model") +
  ylim(0, 40) +
  xlab("Model") +
  ylab("Count of Datasets with Lowest AICc") +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~7)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
aic_resonly_plot_d7


## MRM results ##
extra_arg <- "_multi_unoptimized_scaled"
mrm_results <- read.csv(paste0('mrm_r2_adj_results', extra_arg, '.csv'))

colnames(mrm_results) <- gsub("_dist", "", colnames(mrm_results))
colnames(mrm_results) <- gsub("\\.\\.\\.", "_", colnames(mrm_results))
mrm_resonly_cols_to_keep <- c("biores", "landres", "allres",
                              "geog_biores", "geog_landres", "geog_allres",
                              "env_biores", "env_landres", "env_allres",
                              "geog_env_biores", "geog_env_landres", "geog_env_allres")


## Resistance Only, unoptimized ##
mrm_unoptimized_for_plot <- pivot_longer(mrm_results, 
                                         cols = mrm_resonly_cols_to_keep, 
                                         names_to = 'Model', 
                                         values_to = 'R2')

mrm_unoptimized_for_plot$Model <- factor(mrm_unoptimized_for_plot$Model, levels = mrm_resonly_cols_to_keep)

# Group models by number of predictors for plotting
mrm_unoptimized_for_plot$Group <- sapply(mrm_unoptimized_for_plot$Model, function(x) {
  n <- lengths(regmatches(x, gregexpr("_", x)))
  if (n == 0) {
    "Single"
  } else if (n == 1) {
    "Double"
  } else {
    "Triple"
  }
})
mrm_unoptimized_for_plot$Group <- factor(mrm_unoptimized_for_plot$Group, levels = c("Single", "Double", "Triple"))

#  plot with replaced model labels
mrm_resonly_plot <- ggplot(mrm_unoptimized_for_plot, aes(x = Model, y = R2, fill = Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Reds") +
  xlab("Model") +
  ylab(expression("Adjusted " * italic(R)^2)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
mrm_resonly_plot

colMeans(mrm_results[,mrm_resonly_cols_to_keep])

# MRM rankings
mrm_resonly_rankings <- mrm_results[,mrm_resonly_cols_to_keep]
mrm_resonly_ranks <- t(apply(-mrm_resonly_rankings, 1, rank))

mrm_resonly_ranksums <- colSums(mrm_resonly_ranks ==1)

# Convert to a data frame for ggplot
mrm_df_first_res <- data.frame(
  Model = names(mrm_resonly_ranksums),
  Count = as.numeric(mrm_resonly_ranksums)
)
mrm_df_first_res$Model <- factor(mrm_df_first_res$Model, mrm_resonly_cols_to_keep)

mrm_df_first_res$Group <- sapply(mrm_df_first_res$Model, function(x) {
  n <- lengths(regmatches(x, gregexpr("_", x)))
  if (n == 0) {
    "Single"
  } else if (n == 1) {
    "Double"
  } else {
    "Triple"
  }
})
mrm_df_first_res$Group <- factor(mrm_df_first_res$Group, levels = c("Single", "Double", "Triple"))


# Plot the counts with model labels
mrm_resonly_plot_count <- ggplot(mrm_df_first_res, aes(x = Model, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Reds") +
  xlab("Model") +
  ylab(expression("Count of Datasets with Highest Adjusted " * italic(R)^2)) +
  ggtitle("MRM") +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 8)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
mrm_resonly_plot_count

# find the top ranked driver for each species
mrm_ranks_df <- data.frame(mrm_resonly_ranks)
mrm_ranks_df$top_ranked <- apply(mrm_ranks_df[, mrm_resonly_cols_to_keep], 1, function(x) mrm_resonly_cols_to_keep[which.min(x)])
mrm_ranks_df$top_ranked <- factor(mrm_ranks_df$top_ranked, levels = mrm_resonly_cols_to_keep)
# count number of times each distance is included in best model
mrm_ranks_df$has_biores    <- ifelse(grepl("biores", mrm_ranks_df$top_ranked), 1, 0)
mrm_ranks_df$has_landres    <- ifelse(grepl("landres", mrm_ranks_df$top_ranked), 1, 0)
mrm_ranks_df$has_allres <- ifelse(grepl("allres", mrm_ranks_df$top_ranked), 1, 0)

has_res <- data.frame(count=colSums(mrm_ranks_df[,c('has_biores', 'has_landres', 'has_allres')]))
has_res$Model <- gsub("has_", "", rownames(has_res))
has_res$Model <- factor(has_res$Model, levels = c("biores", "landres", "allres"))
has_res_plot <- ggplot(has_res, aes(x=Model, y=count, fill=Model)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=count), vjust=-0.5, size=3) +
  ylim(c(0,20)) +
  scale_fill_brewer(palette = "Reds") +
  xlab("Model") +
  ylab("Count") +
  ggtitle("MRM") +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
has_res_plot

# PLOT TOGETHER
aic_resonly_plot_d2 <- aic_resonly_plot_d2
aic_resonly_plot_d7 <- aic_resonly_plot_d7 + ylab(NULL)
has_res_plot <- has_res_plot + ylab(NULL)
top_bar <- ggarrange(aic_resonly_plot_d2, aic_resonly_plot_d7,
                     ncol = 2, nrow = 1, 
                     labels = c("A", "B"),
                     label.x = 0.05)
bottom_bar <- ggarrange(mrm_resonly_plot_count, has_res_plot,
                        ncol = 2, nrow = 1, 
                        labels = c("C", "D"),
                        label.x = 0.05,
                        align = "h")
resonly_combined_plot <- ggarrange(top_bar, bottom_bar, 
                                   ncol = 1, heights = c(1, 1.5))

# Display the final plot
print(resonly_combined_plot)

ggsave("plots/resonly_plots.svg", plot = resonly_combined_plot, width = 10, height = 8)

##########################################################
## Optimized Landscape Resistance Only against IBD, IBE ##
##########################################################
extra_arg <- "_multi_optimized_scaled_AICc"
modresults_landonly <- read.csv(paste0('aic_results', extra_arg, '.csv'))

# Keep optimized landscape resistance only
cols_to_keep_landonly <- c("geog_dist", "env_dist", "landres_dist")
mod_order_landonly <- c("IBD", "IBE", "IBR[LAND]")

# Rank models, allowing multiple best ranked models by delta AICc threshold
aic_landonly <- as.matrix(modresults_landonly[,cols_to_keep_landonly])
aic_landonly_adjusted_ranks_d2 <- t(apply(aic_landonly, 1,function(x) adjusted_rank(x, delta_threshold=2)))
colnames(aic_landonly_adjusted_ranks_d2) <- mod_order_landonly
aic_landonly_ranksums_adj_d2 <- colSums(aic_landonly_adjusted_ranks_d2 ==1)
aic_landonly_ranksums_adj_d2

# Convert to a data frame for ggplot
df_first_land_d2 <- data.frame(
  Model = names(aic_landonly_ranksums_adj_d2),
  Count = as.numeric(aic_landonly_ranksums_adj_d2)
)
df_first_land_d2$Model <- factor(df_first_land_d2$Model, mod_order_landonly)

# Plot the counts with model labels mapped
aic_landonly_plot_d2 <- ggplot(df_first_land_d2, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  ylim(0, 40) +
  scale_fill_brewer(palette = "Blues") +
  ylab("Count of Datasets with Lowest AICc") +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~2)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
aic_landonly_plot_d2

## repeat for delta 7 ##
aic_landonly_adjusted_ranks_d7 <- t(apply(aic_landonly, 1,function(x) adjusted_rank(x, delta_threshold=7)))
colnames(aic_landonly_adjusted_ranks_d7) <- mod_order_landonly
aic_landonly_ranksums_adj_d7 <- colSums(aic_landonly_adjusted_ranks_d7 ==1)
aic_landonly_ranksums_adj_d7

# Convert to a data frame for ggplot
df_first_land_d7 <- data.frame(
  Model = names(aic_landonly_ranksums_adj_d7),
  Count = as.numeric(aic_landonly_ranksums_adj_d7)
)
df_first_land_d7$Model <- factor(df_first_land_d7$Model, mod_order_landonly)

# Plot the counts with model labels mapped
aic_landonly_plot_d7 <- ggplot(df_first_land_d7, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  ylim(0, 40) +
  scale_fill_brewer(palette = "Blues") +
  ylab("Count of Datasets with Lowest AICc") +
  ggtitle(expression(LME - Delta~AICc~"\u2264"~7)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
aic_landonly_plot_d7


## MRM ##
extra_arg <- "_multi_optimized_scaled"
mrm_results_landonly <- read.csv(paste0('mrm_r2_adj_results', extra_arg, '.csv'))

colnames(mrm_results_landonly) <- gsub("_dist", "", colnames(mrm_results_landonly))
colnames(mrm_results_landonly) <- gsub("\\.\\.\\.", "_", colnames(mrm_results_landonly))
mrm_landonly_cols_to_keep <- c("geog", "env", "landres",
                               "geog_env", "geog_landres", "env_landres",
                               "geog_env_landres")
mrm_landonly <- mrm_results_landonly[,mrm_landonly_cols_to_keep]
mrm_landonly_ranks <- t(apply(-mrm_landonly, 1, rank))

mrm_landonly_ranksums <- colSums(mrm_landonly_ranks ==1)

# Melting the data for plotting
# Convert to a data frame for ggplot
mrm_df_first_land <- data.frame(
  Model = names(mrm_landonly_ranksums),
  Count = as.numeric(mrm_landonly_ranksums)
)
mrm_df_first_land$Model <- factor(mrm_df_first_land$Model, mrm_landonly_cols_to_keep)

mrm_landonly_plot <- ggplot(mrm_df_first_land, aes(x=Model, y=Count, fill=Model))+
  geom_bar(stat = 'identity', position="stack")+
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette='Reds',name="IBR Model")+
  ylim(c(0,12))+
  xlab("Model")+
  ggtitle("MRM") +
  ylab(expression("Count of Datasets with Highest Adjusted " * italic(R)^2)) +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))
mrm_landonly_plot

mrm_landonly_for_plot <- pivot_longer(mrm_landonly, 
                                      cols = mrm_landonly_cols_to_keep, 
                                      names_to = 'Model', 
                                      values_to = 'R2')

mrm_landonly_for_plot$Model <- factor(mrm_landonly_for_plot$Model, 
                                      levels = mrm_landonly_cols_to_keep)

mrm_landonly_plot_r2 <- ggplot(mrm_landonly_for_plot, aes(x=Model, y=R2, fill=Model))+
  geom_boxplot()+
  scale_fill_brewer(palette='Reds')+
  xlab("Model")+
  ylab(expression("Adjusted " * italic(R)^2))+
  ggtitle("MRM")+
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 12))
mrm_landonly_plot_r2


mrm_ranks_df_landonly <- data.frame(mrm_landonly_ranks)
mrm_ranks_df_landonly$top_ranked <- apply(mrm_ranks_df_landonly[, mrm_landonly_cols_to_keep], 1, function(x) mrm_landonly_cols_to_keep[which.min(x)])
mrm_ranks_df_landonly$top_ranked <- factor(mrm_ranks_df_landonly$top_ranked, levels = mrm_landonly_cols_to_keep)
# count number of times each distance is included in best model
mrm_ranks_df_landonly$has_geog    <- ifelse(grepl("geog", mrm_ranks_df_landonly$top_ranked), 1, 0)
mrm_ranks_df_landonly$has_env    <- ifelse(grepl("env", mrm_ranks_df_landonly$top_ranked), 1, 0)
mrm_ranks_df_landonly$has_landres <- ifelse(grepl("landres", mrm_ranks_df_landonly$top_ranked), 1, 0)

has_mod <- data.frame(count=colSums(mrm_ranks_df_landonly[,c('has_geog', 'has_env', 'has_landres')]))
has_mod$Model <- gsub("has_", "", rownames(has_mod))
has_mod$Model <- factor(has_mod$Model, levels = c("geog", "env", "landres"))

has_mod_plot <- ggplot(has_mod, aes(x=Model, y=count, fill=Model)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=count), vjust=-0.5, size=3) +
  ylim(c(0,40)) +
  scale_fill_brewer(palette = "Reds") +
  xlab("Model") +
  ylab("Count") +
  ggtitle("MRM") +
  scale_x_discrete(labels = function(x) sapply(x, map_model_label)) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")
has_mod_plot

## Arrange graphs
aic_landonly_plot_d7 <- aic_landonly_plot_d7 + ylab(NULL)
has_mod_plot <- has_mod_plot + ylab(NULL)
top_row <- ggarrange(aic_landonly_plot_d2, aic_landonly_plot_d7,
                     ncol=2, widths = c(1, 1, 0.5),
                     labels=c("A", "B"),
                     label.x = 0.05)
bottom_row <- ggarrange(mrm_landonly_plot,has_mod_plot,
                        labels=c("C", "D"),
                        label.x = 0.05,
                        # label.y = 1.05,
                        align="h")
driver_plot <- ggarrange(top_row, bottom_row, 
                         ncol=1, nrow=2)

driver_plot

ggsave("plots/landonly_plots.svg", plot = driver_plot, width = 10, height = 8)


#########################
## Variable Importance ##
#########################
varimp <- read.csv('variable_importance.csv')
varimp_pc <- read.csv('rel_var_importance.csv')

temp <- c("bio1", "bio2", "bio4", "bio8")
precip <- c("bio12", "bio15", "bio18")
human <- c("Croplands", "Pasture", "Built")
veg <- c("Tree", "Shrub", "Grass")

temp_color <- "#E74C3C"
precip_color <- "#4B9CD3"
human_color <- "#7F8C8D"
veg_color <- "#2ECC71"

# Combine colors into a named vector for the legend
legend_colors <- c("Temperature" = temp_color,
                   "Precipitation" = precip_color,
                   "Human" = human_color,
                   "Vegetation" = veg_color)

# Create a mapping from variable to type
variable_to_type <- c(rep("Temperature", length(temp)),
                      rep("Precipitation", length(precip)),
                      rep("Human", length(human)),
                      rep("Vegetation", length(veg)))
names(variable_to_type) <- c(temp, precip, human, veg)

# Calculate median and reorder variables by median
calculate_median_order <- function(df) {
  medians <- apply(df, 2, median, na.rm = TRUE)
  order <- names(sort(medians, decreasing = TRUE))
  return(order)
}

bio_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'bioclim', 9:15])
land_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'landscape', 3:8])
all_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'all', 3:15])

# BIOCLIM model
bio_varimp_pc <- varimp_pc[varimp_pc$model=='bioclim',]
bio_varimp_melted <- melt(bio_varimp_pc[,9:15])
bio_varimp_melted$variable_type <- variable_to_type[as.character(bio_varimp_melted$variable)]
bio_varimp_boxplot <- ggplot(bio_varimp_melted, aes(x=reorder(variable, value, FUN = median), y=value, fill=variable_type))+
  geom_boxplot()+
  scale_fill_manual(values = legend_colors) +
  coord_flip()+
  ylim(c(0,0.9))+
  ylab('Relative Importance')+
  xlab('Variable')+
  theme_minimal()+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
bio_varimp_boxplot
# LANDSCAPE model
land_varimp_pc <- varimp_pc[varimp_pc$model=='landscape',]
land_varimp_melted <- melt(land_varimp_pc[,3:8])
land_varimp_melted$variable_type <- variable_to_type[as.character(land_varimp_melted$variable)]
land_varimp_boxplot <- ggplot(land_varimp_melted, aes(x=reorder(variable, value, FUN = median), y=value, fill=variable_type))+
  geom_boxplot()+
  scale_fill_manual(values = legend_colors) +
  coord_flip()+
  ylim(c(0,0.9))+
  ylab('Relative Importance')+
  xlab('Variable')+
  theme_minimal()+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
land_varimp_boxplot
# ALL model
all_varimp_pc <- varimp_pc[varimp_pc$model=='all',]
all_varimp_melted <- melt(all_varimp_pc[,3:15])
all_varimp_melted$variable_type <- variable_to_type[as.character(all_varimp_melted$variable)]
all_varimp_boxplot <- ggplot(all_varimp_melted, aes(x=reorder(variable, value, FUN = median), y=value, fill=variable_type))+
  geom_boxplot()+
  scale_fill_manual(values = legend_colors) +
  coord_flip()+
  ylim(c(0, 0.9))+
  ylab('Relative Importance')+
  xlab('Variable')+
  theme_minimal()+
  guides(fill = guide_legend(title = "Variable Type")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
all_varimp_boxplot

# Arrange the plots
combined_plot <- plot_grid(
  plot_grid(bio_varimp_boxplot, land_varimp_boxplot, ncol = 2),
  all_varimp_boxplot,
  ncol = 1,
  rel_heights = c(1, 1.5)  # Adjust the relative heights as needed
)

# Display the final plot
print(combined_plot)