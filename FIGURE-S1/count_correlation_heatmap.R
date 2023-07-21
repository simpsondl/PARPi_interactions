library(reshape2)
library(ggplot2)


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# Import data
counts <- read_tsv("../DATA/raw_counts_all_arms.txt")

# Get correlations
corcounts <- round(cor(counts[,c("T0.R1", "T0.R2",
                                 "DMSO.R1", "DMSO.R2",
                                 "NIRAP.R1", "NIRAP.R2")]), 2)
# Remove redundancy
upper.corcounts <- get_upper_tri(corcounts)
# Melt for plotting
corcounts_melt <- melt(upper.corcounts, na.rm = TRUE)

# Make heatmap
corplot <- ggplot(corcounts_melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(0.25,1), space = "Lab",
                      name = "Pearson Correlation") +
  geom_text(aes(Var2, Var1, label = value), color = "black") +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1,0),
        legend.position = c(.4,.7),
        legend.direction = "horizontal",
        plot.margin = margin(l = 1, unit = "cm")) +
  guides(fill = guide_colorbar(barwidth = 7, 
                               barheight = 1,
                               title.position = "top",
                               title.hjust = 0.5)) +
  coord_fixed()

ggsave("../FIGURES/FIGURE-S1/library_distribution_count_correlations_heatmap_noplasmid.svg", corplot, 
       device = "svg", width = 4, height = 4)
