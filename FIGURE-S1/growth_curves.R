library(ggplot2)
library(reshape2)

growth_curves <- read.delim("../../DATA/External_Data/growth_curves.txt")
growth.melt <- melt(growth_curves, id.vars = "X")
growth.melt$X <- gsub("day ", "", growth.melt$X)

plt <- ggplot(growth.melt, aes(X, value, col = variable, group = variable)) + 
  geom_line() + 
  theme_bw() + 
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = "Days Post Dosing", y = "Population Doublings", col = "Arm") +
  scale_color_manual(values = c("#999999", "#777777", "darkorchid4", "darkorchid3"))

ggsave("../../FIGURES/FIGURE-S1/growth_curve.svg", plt, device = "svg", height = 4, width = 5)
