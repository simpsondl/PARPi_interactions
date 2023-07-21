library(readr)
library(dplyr)
library(ggplot2)

gamma.bound.sig <- 1.5274 #qvalue <= .05
tau.bound.sig <- 1.5301 #qvalue <= .05
nu.bound.sig <- 2.977 #qvalue <= .20

gi.cmp2 <- read_tsv("../DATA/Interaction_Scores/InteractionCategories/genecombination_interactioncategories_all_interactions.txt")

gi.cmp2$InteractionCategory <- NA
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & abs(gi.cmp2$Tau) < tau.bound.sig] <- "No interaction"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau > 0] <- "Interaction, unmodified - buffering"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau < 0] <- "Interaction, unmodified - synthetic lethal"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau > 0] <- "Interaction, reversed - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau < 0] <- "Interaction, reversed - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau > 0] <- "Interaction, possibly reversed - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau < 0] <- "Interaction, possibly reversed - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, novel - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, novel - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, possibly novel - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau <= -1*tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, possibly novel - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig] <- "Interaction, masked - gamma positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig] <- "Interaction, masked - gamma negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, possibly masked - gamma positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, possibly masked - gamma negative"
gi.cmp2$InteractionCategory[gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, modified - buffering, nu positive"
gi.cmp2$InteractionCategory[gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, modified - buffering, nu negative"
gi.cmp2$InteractionCategory[gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, modified - synthetic lethal, nu positive"
gi.cmp2$InteractionCategory[gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, modified - synthetic lethal, nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & abs(gi.cmp2$Tau) < tau.bound.sig & abs(gi.cmp2$Nu) > nu.bound.sig] <- "Uncharacterized"

gi.cmp2$InteractionCategory <- factor(gi.cmp2$InteractionCategory, 
                                      levels = unique(gi.cmp2$InteractionCategory)[c(6,
                                                                                     1,5,15,17,
                                                                                     11,10,7,14,
                                                                                     3,2,9,
                                                                                     18,12,19,
                                                                                     4,8,20,16,
                                                                                     13)])

gi.samp <- rbind(gi.cmp2[gi.cmp2$InteractionCategory != "No interaction",],
                 gi.cmp2[gi.cmp2$InteractionCategory == "No interaction",][1:1000,])

plt <- ggplot(gi.cmp2, aes(Gamma, Tau, col = InteractionCategory)) + 
  theme_bw() + 
  scale_color_manual(values = c("#777777",
                                         "#051e39", "#69abf2", "#69abf2", "#051e39",
                                         "darkorchid4", "#c68ce3", "#c68ce3", "darkorchid4",
                                         "#aa2244", "#eeaabb", "#aa2244",
                                         "#aa2244", "#eeaabb", "#aa2244",
                                         "#346f66", "#acd8d1", "#acd8d1", "#346f66",
                                "#777777")) +
                                           geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound.sig), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .8, size = .8) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "None") + coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("") + ylab("")

plt_lab <- ggplot(gi.samp, aes(Gamma, Tau, col = InteractionCategory)) + 
  theme_bw() + 
  scale_color_manual(values = c("#777777",
                                         "#051e39", "#69abf2", "#69abf2", "#051e39",
                                         "darkorchid4", "#c68ce3", "#c68ce3", "darkorchid4",
                                         "#aa2244", "#eeaabb", "#aa2244",
                                         "#aa2244", "#eeaabb", "#aa2244",
                                         "#346f66", "#acd8d1", "#acd8d1", "#346f66",
                                "#777777")) +
                                           geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound.sig, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound.sig), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .8, size = .8) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) + coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("Gamma IS") + ylab("Tau IS") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave("../FIGURES/FIGURE-S5/interaction_category_definitions.png", plt, device = "png",
       height = 10, width = 6, dpi = 300)
ggsave("../FIGURES/FIGURE-S5/interaction_category_definitions_labels.svg", plt_lab, device = "svg",
       height = 10, width = 8, dpi = 300)
