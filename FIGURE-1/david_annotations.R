library(readr)
library(igraph)
library(ggraph)

david2 <- read_tsv("../../DATA/External_Data/david2.txt")

ddr <- unique(c(strsplit(david2$Genes[1], ", ")[[1]], strsplit(david2$Genes[2], ", ")[[1]], strsplit(david2$Genes[3], ", ")[[1]]))

replication <- unique(c(unlist(strsplit(david2$Genes[4:7], ", "))))

cc <- unique(c(unlist(strsplit(david2$Genes[10:13], ", "))))

transport <- unique(c(unlist(strsplit(david2$Genes[15:21], ", "))))

senescense <- unique(c(unlist(strsplit(david2$Genes[22:24], ", "))))

biorhythms <- unique(c(unlist(strsplit(david2$Genes[34:36], ", "))))

mrnaprocessing <- unique(c(unlist(strsplit(david2$Genes[62:67], ", "))))

immune <- unique(c(unlist(strsplit(david2$Genes[c(68:79)], ", "))))

signaling <- unique(c(unlist(strsplit(david2$Genes[215:224], ", "))))

#======================
david3 <- read_tsv("../../DATA/External_Data/david3.txt")

baseexcisionrepair <- unique(c(unlist(strsplit(david3$Genes[c(122,172,181)], ", "))))

nucleotideexcisionrepair <- unique(c(unlist(strsplit(david3$Genes[c(133,158,640)], ", "))))

mismatchrepair <- unique(c(unlist(strsplit(david3$Genes[c(80,895)], ", "))))

nhej <- unique(c(unlist(strsplit(david3$Genes[c(42,156,162,336)], ", "))))

hr <- unique(c(unlist(strsplit(david3$Genes[c(13,18,97,216,297)], ", "))))

interstrandcrosslink <- unique(c(unlist(strsplit(david3$Genes[c(33)], ", "))))

telomere <- unique(c(unlist(strsplit(david3$Genes[grepl("telomere maintenance", david3$Term, ignore.case = TRUE)], ", "))))

###################
# Circle packing
###################

edges <- data.frame(from = "meta",
                    to = c("meta.DNADamageRepair", "meta.DNAReplication", "meta.CellCycle", "meta.Transport", "meta.CellularSenescence",
                           "meta.Biorhythms", "meta.mRNAProcessing", "meta.ImmuneResponse", "meta.Signaling"))
edges <- rbind(edges,
               data.frame(from = "meta.DNADamageRepair",
                          to = c("meta.DNADamageRepair.BER", "meta.DNADamageRepair.NER", "meta.DNADamageRepair.MMR",
                                 "meta.DNADamageRepair.NHEJ", "meta.DNADamageRepair.HR", 
                                 "meta.DNADamageRepair.InterstrandCross-link", "meta.DNADamageRepair.TelomereMaintenance")))
vertices <- data.frame(name = unique(c(edges$from,edges$to)), 
                       size = c(0, 141, 57, 90, 38, 45, 24, 34, 55, 23, 
                                20, 26, 13, 29, 67, 22, 48))
vertices$shortName <- gsub(".*\\.","",vertices$name)
vertices$shortName[vertices$shortName == "meta"] <- NA

mygraph <- graph_from_data_frame( edges, vertices=vertices )

p <- ggraph(mygraph, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = as.factor(depth), col = as.factor(depth))) +
  theme_void() + 
  theme(legend.position="FALSE") +
  geom_node_label(aes(label = shortName), repel = FALSE) +
  scale_fill_manual(values=c("0" = "white", "1" = "#a9a9ad", "2" = "#3083bf")) +
  scale_color_manual(values=c("0" = "white", "1" = "#a9a9ad", "2" = "#3083bf"))

# Organize for saving
annos <- data.frame(Gene = c(ddr, replication, cc, transport,
                             senescense, biorhythms, mrnaprocessing, immune, signaling,
                             baseexcisionrepair, nucleotideexcisionrepair, mismatchrepair, nhej, 
                             hr, interstrandcrosslink, telomere), 
                    Annotation = c(rep("DNA Damage and Repair", length(ddr)),
                                   rep("DNA Replication", length(replication)),
                                   rep("Cell Cycle", length(cc)),
                                   rep("mRNA Transport", length(transport)),
                                   rep("Cellular Senescence", length(senescense)),
                                   rep("Biorhythms", length(biorhythms)),
                                   rep("mRNA Processing", length(mrnaprocessing)),
                                   rep("Immune", length(immune)),
                                   rep("Signaling", length(signaling)),
                                   rep("Base Excision Repair", length(baseexcisionrepair)),
                                   rep("Nucleotide Excision Repair", length(nucleotideexcisionrepair)),
                                   rep("Mismatch Repair", length(mismatchrepair)),
                                   rep("Nonhomologous end joining", length(nhej)),
                                   rep("Homologous Recombination", length(hr)),
                                   rep("Interstrand Cross-link Repair", length(interstrandcrosslink)),
                                   rep("Telomere Maintenance", length(telomere))))

ggsave("../FIGURES/FIGURE-1/library_annotation_composition.svg", p, device = "svg", width = 6, height = 6, dpi = 300)
write_tsv(annos, "../../DATA/External_Data/david_assigned_annotations.txt")
