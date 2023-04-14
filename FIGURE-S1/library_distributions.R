library(readr)
library(ineq)
library(reshape2)
library(wesanderson)
library(ggplot2)

counts <- read_tsv("../../DATA/raw_counts_all_arms.txt")

######## DEFINE DATA FOR LORENZ CURVES ###########
lc.t0 <- Lc(counts$T0.R1)
lc.dmso <- Lc(counts$DMSO.R1)
lc.nirap <- Lc(counts$NIRAP.R1)
lc.plasmid <- Lc(counts$PLASMID)
tmp.df <- data.frame("Percent.Data" = lc.t0$p, 
                     "T0" = lc.t0$L,
                     "NIRAP" = lc.nirap$L,
                     "PLASMID" = lc.plasmid$L)

tmp.df2 <- tmp.df[seq(1,1311026,131),]
tmp.melt <- melt(tmp.df2, id.vars = "Percent.Data")
colnames(tmp.melt)[2] <- "Library"

lz <- ggplot(tmp.melt, aes(Percent.Data, value, color = Library)) + 
  scale_color_manual(values = c("#44AA99", "darkorchid4", "black")) +
  geom_line() + 
  theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "#cccccc") +
  xlab("Percent of Library Elements") +
  ylab("Percent of Reads") +
  ggtitle("Library Skew") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  annotate("text", x = .05, y = .9, hjust = 0, label = "Gini coefficients:") +
  annotate("text", hjust = 0, x= .1, y = .82, label = paste("PLASMID =", round(Gini(counts$PLASMID),3))) +
  annotate("text", hjust = 0, x= .1, y = .76, label = paste("T0 =", round(Gini(counts$T0.R1),3))) +
  annotate("text", hjust = 0, x= .1, y = .7, label = paste("NIRAP =", round(Gini(counts$NIRAP.R1),3))) 

######## DEFINE DATA FOR DENSITY CURVES ###########
counts.norm <- counts[,c(10,14,16)]
counts.norm$T0.R1 <- counts.norm$T0.R1 * sum(counts.norm$NIRAP.R1)/sum(counts.norm$T0.R1)
counts.norm$PLASMID <- counts.norm$PLASMID * sum(counts.norm$NIRAP.R1)/sum(counts.norm$PLASMID)
counts.norm$T0.R1 <- log2(counts.norm$T0.R1)
counts.norm$NIRAP.R1 <- log2(counts.norm$NIRAP.R1)
counts.norm$PLASMID <- log2(counts.norm$PLASMID)
counts.melt <- melt(counts.norm)

plasmid.quant <- quantile(counts.norm$PLASMID, c(.05, .95))
t0.quant <- quantile(counts.norm$T0.R1, c(.05, .95))
nirap.quant <- quantile(counts.norm$NIRAP.R1, c(.05, .95))

dist.dens <- ggplot(counts.melt, aes(value, col = variable)) + 
  scale_color_manual(values = c("#44AA99", "darkorchid4", "black")) +
  geom_density(size = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = expression("log"["2"] ~ "(counts)"),
       y = "density",
       col = "Library") +
  annotate("text", x = 1, hjust = 0, y = .45, label = "95/5 Ratio:") +
  annotate("text", x = 1.2, hjust = 0, y = .40, label = paste("PLASMID =", round(2^(plasmid.quant[[2]] - plasmid.quant[[1]]), 3))) +
  annotate("text", x = 1.2, hjust = 0, y = .37, label = paste("T0 =", round(2^(t0.quant[[2]] - t0.quant[[1]]), 3))) +
  annotate("text", x = 1.2, hjust = 0, y = .34, label = paste("NIRAP =", round(2^(nirap.quant[[2]] - nirap.quant[[1]]), 3)))   

######## DEFINE DATA FOR CORRELATION HEATMAP ###########
corcounts <- round(cor(counts[,10:16]), 2)
# Remove redundancy
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper.corcounts <- get_upper_tri(corcounts)
# Melt for plotting
corcounts_melt <- melt(upper.corcounts, na.rm = TRUE)
corcounts_melt2 <- corcounts_melt[corcounts_melt$Var2 != "PLASMID",]

# Make heatmap
corplot <- ggplot(corcounts_melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(0,1), space = "Lab",
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

# Make heatmap
corplot2 <- ggplot(corcounts_melt2, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient(low = "white", high = "red",
                      limit = c(0,1), space = "Lab",
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


######## SAVE PLOTS ###########
ggsave("../../FIGURES/FIGURE-S1/library_distribution_lorenz_curves.svg", lz, 
       device = "svg", width = 5, height = 4)
ggsave("../../FIGURES/FIGURE-S1/library_distribution_density_curves.svg", dist.dens, 
       device = "svg", width = 5, height = 4)
ggsave("../../FIGURES/FIGURE-S1/library_distribution_count_correlations_heatmap.svg", corplot, 
       device = "svg", width = 4, height = 4)
ggsave("../../FIGURES/FIGURE-S1/library_distribution_count_correlations_heatmap_noplasmid.svg", corplot2, 
       device = "svg", width = 4, height = 4)

