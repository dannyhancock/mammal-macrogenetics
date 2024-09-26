library(ggplot2)
library(readxl)
library(plyr)
library(ggpubr)
library(data.table)
library(gridExtra)

#########################
### MODELLING RESULTS ###
#########################

## AIC results ##
modresults <- read.csv('aic_results.csv')

## resistance only ##
aic_resonly <- modresults[,c('biores','landres','allres')]
aic_resonly_ranks <- t(apply(aic_resonly, 1, rank))
aic_resonly_ranksums <- colSums(aic_resonly_ranks ==1)

# Creating contingency table
aic_cont_resonly <- table(c(col(aic_resonly_ranks)), c(aic_resonly_ranks))
rownames(aic_cont_resonly) <- colnames(aic_resonly_ranks)

# Melting the data for plotting
aic_resonly_for_plot <- melt(aic_cont_resonly)
aic_resonly_for_plot$Var1 <- factor(aic_resonly_for_plot$Var1, 
                                    levels=c('allres', 'biores', 'landres'),
                                    labels=c('ALL', 'BIO', 'LAND'))

aic_resonly_plot <- ggplot(aic_resonly_for_plot, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat = 'identity', position="stack")+
  scale_fill_brewer(palette='Blues',name="IBR Model")+
  ylim(c(0,41))+
  xlab("Rank")+
  ylab("Number of Species")+
  theme_minimal()+
  ggtitle('LME')+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12),
        axis.title = element_text(size=14, face="bold"),
        title = element_text(size=14, face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
aic_resonly_plot

## LAND ONLY ##
# Read and process the data
aic_landonly <- modresults[, c('landres', 'env', 'geog')]
aic_landonly_ranks <- t(apply(aic_landonly, 1, rank))
aic_landonly_ranksums <- colSums(aic_landonly_ranks == 1)

# Creating contingency table
aic_cont_landonly <- table(c(col(aic_landonly_ranks)), c(aic_landonly_ranks))
rownames(aic_cont_landonly) <- colnames(aic_landonly_ranks)

# Melting the data for plotting
aic_landonly_for_plot <- melt(aic_cont_landonly)

# Rename factors directly during assignment to ensure proper order
aic_landonly_for_plot$Var1 <- factor(aic_landonly_for_plot$Var1, 
                                     levels = c("geog", "env", "landres"),
                                     labels = c("IBD", "IBE", "IBR(LAND)"))

# Plotting
aic_landonly_plot <- ggplot(aic_landonly_for_plot, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = 'identity', position = "stack") +
  scale_fill_brewer(palette = 'Reds', name = "Driver") +
  ylim(c(0, 41)) +
  xlab("Rank of Driver") +
  ylab("Number of Species") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
aic_landonly_plot


## MRM results ##
mrm_results <- read.csv('mrm_r2_results.csv')

## Resistance Only ##
# Read and process the data
mrm_resonly <- mrm_results[,c('biores','landres','allres')]
mrm_resonly_ranks <- t(apply(-mrm_resonly, 1, rank))
mrm_resonly_ranksums <- colSums(mrm_resonly_ranks ==1)

# Creating contingency table
mrm_resonly_cont <- table(c(col(mrm_resonly_ranks)), c(mrm_resonly_ranks))
rownames(mrm_resonly_cont) <- colnames(mrm_resonly_ranks)

levels(aic_resonly_for_plot$Var1) <- c("IBRALL", "IBRBIO", "IBRLAND")

# Melting the data for plotting
mrm_resonly_for_plot <- melt(mrm_resonly_cont)
mrm_resonly_for_plot$Var1 <- factor(mrm_resonly_for_plot$Var1, 
                                    levels=c('allres', 'biores', 'landres'),
                                    labels=c("ALL", "BIO", "LAND"))


mrm_resonly_plot <- ggplot(mrm_resonly_for_plot, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat = 'identity', position="stack")+
  scale_fill_brewer(palette='Blues',name="IBR Model")+
  ylim(c(0,41))+
  # xlab("Rank")+
  # ylab("Count")+
  theme_minimal()+
  ggtitle("MRM")+
  theme(text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
mrm_resonly_plot


## Land Only ##
mrm_landonly <- mrm_results[,c('landres','env','geog')]
mrm_landonly_ranks <- t(apply(-mrm_landonly, 1, rank))
mrm_landonly_ranksums <- colSums(mrm_landonly_ranks ==1)

mrm_landonly_ranks
mrm_landonly_ranksums
# Create contingency table
mrm_landonly_cont <- table(c(col(mrm_landonly_ranks)), c(mrm_landonly_ranks))
rownames(mrm_landonly_cont) <- colnames(mrm_landonly_ranks)

# Melting the data for plotting
mrm_landonly_for_plot <- melt(mrm_landonly_cont)
mrm_landonly_for_plot$Var1 <- factor(mrm_landonly_for_plot$Var1, 
                                     levels = c("geog", "env", "landres"),
                                     labels = c("IBD", "IBE", "IBR(LAND)"))

mrm_landonly_plot <- ggplot(mrm_landonly_for_plot, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat = 'identity', position="stack")+
  scale_fill_brewer(palette='Reds',name="Driver")+
  ylim(c(0,41))+
  xlab("Rank")+
  ylab("Count")+
  theme_minimal()+
  theme(text = element_text(size = 12),
        axis.title.y = element_blank(),
        title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
mrm_landonly_plot


## Arrange quad graphs
legend_top <- get_legend(mrm_resonly_plot + theme(legend.position="right"))
legend_bottom <- get_legend(mrm_landonly_plot + theme(legend.position="right"))

# Remove legends, gridlines and x- and y-labels from individual plots
aic_landonly_plot <- aic_landonly_plot + theme(legend.position="none",
                                               axis.title.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank())
aic_resonly_plot <- aic_resonly_plot + theme(legend.position="none",
                                             axis.title.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank())
mrm_resonly_plot <- mrm_resonly_plot + theme(legend.position="none",
                                             axis.title.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank())
mrm_landonly_plot <- mrm_landonly_plot + theme(legend.position="none",
                                               axis.title.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank())


# Arrange in 2x2 grid
mrm_resonly_plot 

top_row <- ggarrange(mrm_resonly_plot, aic_resonly_plot,
                     legend_top, ncol=3, widths = c(1, 1, 0.5),
                     labels=c("A", "B"))
bottom_row <- ggarrange(mrm_landonly_plot + theme(plot.margin=margin(15,0,0,0)), 
                        aic_landonly_plot + theme(plot.margin=margin(15,0,0,0)),
                        legend_bottom, ncol=3, widths = c(1, 1, 0.5),
                        labels=c("C", "D"))
driver_plot <- ggarrange(top_row, bottom_row, 
                         ncol=1, nrow=2)

driver_plot <- annotate_figure(driver_plot,
                               bottom = text_grob("Rank",size=14,face="bold", hjust=2.5),
                               left = text_grob("Number of species", size=14, face="bold", rot=90)
)

driver_plot

ggsave("z_submission/plots/driver_rankings.svg", plot = driver_plot, width = 10, height = 8)

#########################
## Variable Importance ##
#########################
varimp <- read.csv('variable_importance.csv')
varimp_pc <- read.csv('rel_var_importance_2.csv')

temp <- c("bio1", "bio2", "bio4", "bio8")
precip <- c("bio12", "bio15", "bio18")
human <- c("Roads", "Croplands", "Pasture", "Railways", "Built")
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

bio_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'bioclim', 3:9])
land_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'landscape', 10:17])
all_order <- calculate_median_order(varimp_pc[varimp_pc$model == 'all', 3:17])

# BIOCLIM model
bio_varimp_pc <- varimp_pc[varimp_pc$model=='bioclim',]
bio_varimp_melted <- melt(bio_varimp_pc[,3:9])
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

# LANDSCAPE model
land_varimp_pc <- varimp_pc[varimp_pc$model=='landscape',]
land_varimp_melted <- melt(land_varimp_pc[,10:17])
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

# ALL model
all_varimp_pc <- varimp_pc[varimp_pc$model=='all',]
all_varimp_melted <- melt(all_varimp_pc[,3:17])
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

# Arrange the plots
combined_plot <- plot_grid(
  plot_grid(bio_varimp_boxplot, land_varimp_boxplot, ncol = 2),
  all_varimp_boxplot,
  ncol = 1,
  rel_heights = c(1, 1.5)  # Adjust the relative heights as needed
)

# Display the final plot
print(combined_plot)