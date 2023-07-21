kphenos <- read_tsv("../DATA/AUNIP_targeted/HTS_JL021/aunip_targeted_library_all_phenotypes.txt")
rphenos <- read_tsv("../DATA/AUNIP_targeted/HTS_JL020/aunip_targeted_library_all_phenotypes.txt")

# Remove redundancy
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Compress information, only want gamma and tau
kphenos2 <- kphenos[,c(1:11, 15:17, 12:14, 21:23,27:29)]
rphenos2 <- rphenos[,c(1:11, 15:17, 12:14, 21:23,27:29)]

# Get correlations
kcorcounts <- round(cor(kphenos2[,9:ncol(kphenos2)]), 2)
rcorcounts <- round(cor(rphenos2[,9:ncol(rphenos2)]), 2)

comb <- kcorcounts
comb[upper.tri(comb)] <- rcorcounts[upper.tri(rcorcounts)]

# Melt for plotting
corcounts_melt <- melt(comb, na.rm = TRUE)
corcounts_melt$value[corcounts_melt$Var1 == corcounts_melt$Var2] <- NA

# Make heatmap
plt <- ggplot(corcounts_melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") + 
  theme_bw() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       limit = c(0,1), space = "Lab",
                       name = "Pearson Correlation") +
  geom_text(aes(Var2, Var1, label = value), 
            color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[1]) - .5,
                xmax = as.numeric(Var1[3]) + .5,
                ymin = as.numeric(Var1[1]) - .5,
                ymax = as.numeric(Var1[3]) + .5),
            fill = NA, color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[4]) - .5,
                xmax = as.numeric(Var1[6]) + .5,
                ymin = as.numeric(Var1[4]) - .5,
                ymax = as.numeric(Var1[6]) + .5),
            fill = NA, color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[7]) - .5,
                xmax = as.numeric(Var1[9]) + .5,
                ymin = as.numeric(Var1[7]) - .5,
                ymax = as.numeric(Var1[9]) + .5),
            fill = NA, color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[10]) - .5,
                xmax = as.numeric(Var1[12]) + .5,
                ymin = as.numeric(Var1[10]) - .5,
                ymax = as.numeric(Var1[12]) + .5),
            fill = NA, color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[13]) - .5,
                xmax = as.numeric(Var1[15]) + .5,
                ymin = as.numeric(Var1[13]) - .5,
                ymax = as.numeric(Var1[15]) + .5),
            fill = NA, color = "black", size = 2) +
  geom_rect(aes(xmin = as.numeric(Var1[4]) - .5,
                xmax = as.numeric(Var1[12]) + .5,
                ymin = as.numeric(Var1[4]) - .5,
                ymax = as.numeric(Var1[12]) + .5),
            fill = NA, color = "black", size = 1.25) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(l = 1, unit = "cm"), 
        legend.position = "none") +
  coord_fixed()

ggsave("../FIGURES/FIGURE-6/targetscreen_correlation_heatmap.svg",
       plt, device = "svg", height = 6, width = 6)

ggsave("../FIGURES/FIGURE-6/targetscreen_correlation_heatmap.pdf",
       plt, device = "pdf", height = 6, width = 6)
