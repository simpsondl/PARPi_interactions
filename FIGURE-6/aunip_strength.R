library(readr)
library(ggplot2)
library(ggridges)

scores <- read_tsv("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/gene_combination_interaction_scores_nu_oi_avg.txt")
idmap <- read_tsv("../DATA/genecombination_id_map.txt")
scores <- inner_join(scores, idmap)
scores <- scores[scores$Gene1 != scores$Gene2,]
scores <- scores[!grepl("non-targeting", scores$GeneCombinationName),]

scores2 <- scores
scores2$Category <- "All"
scores2$Enrichment <- NA

tmp <- scores[grepl("^AUNIP:|:AUNIP$", scores$GeneCombinationName),]
tmp$Category <- "AUNIP"
for(i in 1:nrow(tmp)){tmp$Enrichment[i] <- sum(tmp$InteractionScore <= tmp$InteractionScore[i])/sum(scores$InteractionScore <= tmp$InteractionScore[i])}
scores2 <- rbind(scores2,
                 tmp)

tmp <- scores[grepl("^FANCA:|:FANCA$", scores$GeneCombinationName),]
tmp$Category <- "FANCA"
for(i in 1:nrow(tmp)){tmp$Enrichment[i] <- sum(tmp$InteractionScore <= tmp$InteractionScore[i])/sum(scores$InteractionScore <= tmp$InteractionScore[i])}
scores2 <- rbind(scores2,
                 tmp)

tmp <- scores[grepl("^RAD54L:|:RAD54L$", scores$GeneCombinationName),]
tmp$Category <- "RAD54L"
for(i in 1:nrow(tmp)){tmp$Enrichment[i] <- sum(tmp$InteractionScore <= tmp$InteractionScore[i])/sum(scores$InteractionScore <= tmp$InteractionScore[i])}
scores2 <- rbind(scores2,
                 tmp)

scores2$Category <- factor(scores2$Category, levels = c("AUNIP", "FANCA", "RAD54L", "All"))

plt.info <- ggplot(scores2, 
                   aes(x = InteractionScore, 
                       y = Category)) + 
  geom_density_ridges_gradient(scale = 3, 
                               rel_min_height = .005,
                               jittered_points = TRUE) +
  scale_fill_viridis_c(name = "% Strong interactions", option = "C") +
  theme_ridges(center_axis_labels = TRUE)

enrichments <- unique(ggplot_build(plt.info)$data[[1]])
enrichments <- enrichments[,c(2,1,5,8)]
enrichments$Enrichment <- NA

tmp <- scores[grepl("^AUNIP:|:AUNIP$", scores$GeneCombinationName),]
tmp$Category <- "AUNIP"
for(i in 1:table(enrichments$group)[[2]]){
  enrichments$Enrichment[enrichments$group == 1][i] <- sum(tmp$InteractionScore <= enrichments$x[enrichments$group == 1][i])/sum(scores$InteractionScore <= enrichments$x[enrichments$group == 1][i]) 
}

tmp <- scores[grepl("^FANCA:|:FANCA$", scores$GeneCombinationName),]
tmp$Category <- "FANCA"
for(i in 1:table(enrichments$group)[[2]]){
  enrichments$Enrichment[enrichments$group == 2][i] <- sum(tmp$InteractionScore <= enrichments$x[enrichments$group == 2][i])/sum(scores$InteractionScore <= enrichments$x[enrichments$group == 2][i]) 
}

tmp <- scores[grepl("^RAD54L:|:RAD54L$", scores$GeneCombinationName),]
tmp$Category <- "RAD54L"
for(i in 1:table(enrichments$group)[[2]]){
  enrichments$Enrichment[enrichments$group == 3][i] <- sum(tmp$InteractionScore <= enrichments$x[enrichments$group == 3][i])/sum(scores$InteractionScore <= enrichments$x[enrichments$group == 3][i]) 
}

enrichments$Enrichment[enrichments$Enrichment == "NaN"] <- 0

# plt <- ggplot(scores2[scores2$Category != "All",], 
#               aes(x = Nu, 
#                   y = Category, 
#                   fill = enrichments$Enrichment[enrichments$x == after_stat(x) & enrichments$group == after_stat(y)])) +
#   geom_density_ridges_gradient(scale = 5, 
#                                rel_min_height = .005,
#                                jittered_points = TRUE,
#                                position = position_points_jitter(width = .05, height = 0),
#                                point_shape = "|", point_size = 3, point_alpha = .7, alpha = .7) +
#   scale_fill_viridis_c(name = "% Strong interactions", option = "C") +
#   theme_ridges(center_axis_labels = TRUE)



plt <- ggplot(scores2, 
              aes(x = InteractionScore, 
                  y = Category,
                  fill = enrichments$Enrichment[enrichments$x == after_stat(x) & enrichments$group == after_stat(y)])) + 
  geom_vline(xintercept = -5.15, linetype = "dashed") +
  geom_density_ridges_gradient(scale = 3, 
                               rel_min_height = .005,
                               jittered_points = TRUE,
                               position = position_points_jitter(height = 0),
                               point_shape = "|", point_size = 1.25, point_alpha = .7, alpha = .7) +
  scico::scale_fill_scico(name = "% strong", palette = "lajolla") +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "none") + xlab("") + ylab("") + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

ggsave("../FIGURES/FIGURE-6/aunip_strength_distribution.png", plt, device = "png",
       height = 4, width = 4, dpi = 300)
