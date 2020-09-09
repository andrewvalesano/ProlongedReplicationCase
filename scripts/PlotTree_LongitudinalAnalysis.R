
### Author: Andrew Valesano
### Purpose: Plot the phylogenetic tree with genomes from perisistent positive case.

# ============================= Modules and data ==========================

library(tidyverse)
library(ggtree)
library(wesanderson)

tree <- read.tree("data/raw/context.iqtree.nwk")
pangolin <- read.csv("data/raw/pangolin_results.csv", stringsAsFactors = FALSE) %>%
  mutate(name = gsub("/", "_", Sequence.name)) %>%
  mutate(Lineage = ifelse(str_detect(Sequence.name, "day"), "Patient Sample (Lineage B.1)", Lineage)) %>%
  select(name, Lineage)

# ========================== Plot the tree! ======================

tree.plot <- ggtree(tree, layout = "circular")

palette_D1 <- wes_palette("Darjeeling1")
palette_D2 <- wes_palette("Darjeeling2")
colors <- c("goldenrod1", palette_D1[2], "plum2", "black", palette_D2[2], palette_D1[5], "mediumpurple", palette_D1[1], "white")

tree.plot.meta <- tree.plot %<+% pangolin +
  geom_point(aes(color = Lineage, alpha = ifelse(is.na(Lineage), 0, 1)), size = 2) +
  theme(legend.position = "right") +
  scale_color_manual(values = colors) +
  geom_treescale() +
  geom_tiplab(size = 2, offset = 0.00001) # 10 by 10

labels <- tree.plot$data
labels <- labels[!labels$isTip,]
labels$label <- as.numeric(labels$label)
labels <- labels[labels$label > 70,]

tree.plot.meta.bootstrap <- tree.plot.meta +
  geom_text(data = labels, aes(label = label), size = 3, nudge_x = -0.000003, nudge_y = 1, color = "black")

