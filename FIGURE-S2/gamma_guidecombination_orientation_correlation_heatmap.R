library(readr)
library(ggplot2)
library(reshape2)

# Load data
phenotypes <- read_tsv("../DATA/raw_phenotypes.txt")
rho <- read_tsv("../DATA/raw_rho_phenotypes.txt")
phenotypes <- inner_join(phenotypes,
                         rho,
                         by = colnames(phenotypes)[1:11])
# Rearrange data for plotting
tmp.df <- rbind(phenotypes[phenotypes$Category == "X+Y",], 
                phenotypes[phenotypes$Category == "X+NT",],
                phenotypes[phenotypes$Category == "NT+NT",],
                phenotypes[phenotypes$Category == "X+X",])

# Identify guide combinations with both AB and BA orientations
gc.freq <- as.data.frame(table(tmp.df$GuideCombinationID))
# Rearrange replicate-averaged gamma orientation-specific phenotypes for plotting
# Keep only guide combinations with both orientations available
all.orientation <- pivot_wider(tmp.df[tmp.df$GuideCombinationID %in% gc.freq$Var1[gc.freq$Freq == 2],
                                        c("Category", "Orientation", "GuideCombinationID", 
                                          "Gamma.R1", "Gamma.R2", "Gamma.Avg",
                                          "Tau.R1", "Tau.R2", "Tau.Avg",
                                          "Rho.R1", "Rho.R2", "Rho.Avg")],
                                 names_from = Orientation,
                                 values_from = c(Gamma.R1, Gamma.R2, Gamma.Avg,
                                                 Tau.R1, Tau.R2, Tau.Avg,
                                                 Rho.R1, Rho.R2, Rho.Avg))

g <- all.orientation[,3:8]
t <- all.orientation[,9:14]
r <- all.orientation[,15:20]
# Update column names for plotting
colnames(g) <- c("R1 - AB", "R1 - BA",
                 "R2 - AB", "R2 - BA",
                 "Avg - AB", "Avg - BA")
colnames(t) <- c("R1 - AB", "R1 - BA",
                 "R2 - AB", "R2 - BA",
                 "Avg - AB", "Avg - BA")
colnames(r) <- c("R1 - AB", "R1 - BA",
                 "R2 - AB", "R2 - BA",
                 "Avg - AB", "Avg - BA")

# Get correlations
corcountsg <- round(cor(g), 2)
corcountst <- round(cor(t), 2)
corcountsr <- round(cor(r), 2)

# Remove redundancy
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper.corcountsg <- get_upper_tri(corcountsg)
upper.corcountst <- get_upper_tri(corcountst)
upper.corcountsr <- get_upper_tri(corcountsr)
# Melt for plotting
corcounts_meltg <- melt(upper.corcountsg, na.rm = TRUE)
corcounts_meltt <- melt(upper.corcountst, na.rm = TRUE)
corcounts_meltr <- melt(upper.corcountsr, na.rm = TRUE)

# Make heatmap
gplot <- ggplot(corcounts_meltg, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(.25,1), space = "Lab",
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

tplot <- ggplot(corcounts_meltt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(.25,1), space = "Lab",
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

rplot <- ggplot(corcounts_meltr, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(.25,1), space = "Lab",
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

ggsave("../FIGURES/FIGURE-S2/gamma_orientation_correlation_heatmap.svg", 
       gplot, device = "svg", width = 6, height = 6)
ggsave("../FIGURES/FIGURE-S2/tau_orientation_correlation_heatmap.svg", 
       tplot, device = "svg", width = 6, height = 6)
ggsave("../FIGURES/FIGURE-S2/rho_orientation_correlation_heatmap.svg", 
       rplot, device = "svg", width = 6, height = 6)
