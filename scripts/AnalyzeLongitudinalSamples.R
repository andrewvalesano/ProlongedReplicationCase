
### Author: Andrew Valesano
### Purpose: Analyze longitudinal samples from persistent positive case.

# ========================= Modules and data ============================

library(tidyverse)
library(patchwork)
library(lubridate)

metadata <- read.csv("data/metadata/metadata_final.csv", stringsAsFactors = FALSE)
variants <- read.csv("data/raw/all.variants.csv", stringsAsFactors = FALSE) %>% mutate(mutation = paste0(REF, POS, ALT))
coverage <- read.csv("data/raw/coverage.csv", stringsAsFactors = FALSE) %>% mutate(ID = gsub("/", "_", ID))

# ============================== Examine consensus mutations =======================

variants_fixed <- filter(variants, ALT_FREQ > 0.5)
variants_fixed_meta <- left_join(variants_fixed, select(metadata, ID, day), by = "ID")

# ======================== Coverage ====================

cov.plot <- ggplot(coverage, aes(x = pos, y = cov)) +
  geom_line() +
  facet_wrap(~ID) +
  scale_y_log10()

# ================================= iSNV in multiple samples ================================

# 

# Get iSNV (2% - 90%) that are found in more than one sample.
isnv <- filter(variants, ALT_FREQ > 0.02 & ALT_FREQ < 0.90)
isnv_meta <- left_join(isnv, select(metadata, ID, day), by = "ID")

isnv_meta %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-")) %>%
  filter(ALT_QUAL > 25 & ALT_DP >= 5 & ALT_RV != 0 & ALT_DP != ALT_RV) %>% 
  group_by(mutation) %>%
  filter(n() > 1) %>%
  ungroup() -> isnv_meta_multiple

# For these iSNV, create dataframe including samples that lack those iSNV.
coverage_meta <- left_join(coverage, select(metadata, ID, day), by = "ID")

isnv_meta_multiple_full <- data.frame()
for(mut in unique(isnv_meta_multiple$mutation))
{
  mut_data <- data.frame(day = unique(metadata$day))
  isnv_meta_multiple_mut <- filter(isnv_meta_multiple, mutation == mut) %>% select(ALT_FREQ, day)
  mut_data <- left_join(mut_data, isnv_meta_multiple_mut, by = "day") %>% mutate(mutation = mut)
  
  position <- unique(filter(isnv_meta_multiple, mutation == mut)$POS)
  
  coverage_meta_mut <- filter(coverage_meta, pos == position) %>% filter(cov >= 100)
  
  mut_data_NA <- filter(mut_data, is.na(ALT_FREQ))
  mut_data_NA <- mutate(mut_data_NA, ALT_FREQ = ifelse(day %in% coverage_meta_mut$day, 0, NA))
  
  mut_data_end <- rbind(mut_data_NA, filter(mut_data, !is.na(ALT_FREQ)))
  isnv_meta_multiple_full <- rbind(isnv_meta_multiple_full, mut_data_end)
}

# ========================================== Get nucleotide diversity per sample ==================================

get_pairwise_distance <- function(df)
{
  if(nrow(df) > 1)
  {
    print(df)
    stopifnot(nrow(df) == 1)
  }
  
  coverage <- df$TOTAL_DP
  cov_factor <- (coverage * (coverage - 1))
  
  cov_minor <- df$ALT_DP
  cov_major <- df$REF_DP
  
  minor_factor <- cov_minor * (cov_minor - 1)
  major_factor <- cov_major * (cov_major - 1)
  
  D <- (cov_factor - (major_factor + minor_factor)) / cov_factor
  
  return(data.frame(site = df$POS, D = D))
}

get_nt_diversity <- function(df)
{
  stopifnot(length(unique(df$ID)) == 1) # should be doing one sample at a time
  
  ID <- as.character(df[1,]$ID)
  
  df %>%
    group_by(POS) %>%
    do(get_pairwise_distance(.)) -> distance_sites
  
  pi <- sum(distance_sites$D) / 29903
  
  df_return <- data.frame(ID = ID, pi = pi)
  
  return(df_return)
}

# Get pi per sample
minor <- filter(variants, 
                  !ID %in% c("day29nps") &
                  ALT_FREQ > 0.02 & 
                  ALT_FREQ < 0.5 & 
                  !str_detect(string = ALT, pattern = "\\+") & 
                  !str_detect(string = ALT, pattern = "\\-") & 
                  ALT_DP >= 10 & 
                  ALT_RV != 0 & 
                  ALT_DP != ALT_RV)
minor %>%
  group_by(ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

# Express pi over time 
pi_meta <- left_join(select(metadata, ID, day), nt_diversity, by = "ID")
pi_meta$pi[is.na(pi_meta$pi)] <- 0
pi_meta <- filter(pi_meta, !ID %in% c("day29nps")) %>% mutate(mutation = "Nucleotide Diversity (pi)")
pi_meta <- mutate(pi_meta, pi_adj = pi*7.3e3) # for plotting on dual axes

# Plot
pi.by.day <- ggplot(pi_meta, aes(x = day, y = pi)) +
  geom_point(size = 2.3) +
  theme_bw() +
  xlab("Day of Illness") +
  ylab("Nucleotide Diversity (pi)") +
  geom_line(size = 1) +
  xlim(c(0, 38))

plot_colors <- c(RColorBrewer::brewer.pal(7, name = "Paired"), "black")
isnv_meta_multiple_full$mutation <- factor(isnv_meta_multiple_full$mutation, levels = c("C5183T", "C5178T", "C6033T", "C10702T", "G25644C", "G28239T", "C28253T"))
isnv.within.hosts <- ggplot(isnv_meta_multiple_full, aes(x = day, y = ALT_FREQ, color = mutation)) +
  geom_point(size = 2.3) +
  geom_line(size = 1) +
  xlab("Days") +
  ylab("Frequency") +
  theme_bw() +
  xlim(c(0, max(metadata$day))) +
  scale_color_manual(name = "Mutation", values = plot_colors)

isnv.within.hosts | pi.by.day

# Plot with dual y-axes!

plot_colors <- c(RColorBrewer::brewer.pal(7, name = "Paired"), "black")
isnv_meta_multiple_full$mutation <- factor(isnv_meta_multiple_full$mutation, levels = c("C5183T", "C5178T", "C6033T", "C10702T", "G25644C", "G28239T", "C28253T"))

isnv.within.hosts <- ggplot(isnv_meta_multiple_full, aes(x = day, y = ALT_FREQ, color = mutation)) +
  geom_point(size = 2.3) +
  geom_line(size = 1) +
  xlab("Day of Illness") +
  ylab("Frequency") +
  theme_bw() +
  xlim(c(0, 40)) +
  scale_color_manual(name = "Mutation", values = plot_colors) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), sec.axis = sec_axis(trans = ~ ./7.3e3, breaks = seq(0, 8e-05, 2e-05), name = "Nucleotide Diversity (pi)")) +
  geom_line(data = pi_meta, aes(x = day, y = pi_adj), size = 1, color = "black") +
  geom_point(data = pi_meta, aes(x = day, y = pi_adj), size = 2.3, color = "black")
isnv.within.hosts # 4 by 6

# ==================================== Plot consensus substitutions ==============================

consensus_data <- read.table("data/raw/samples.nextclade.csv")
header <- as.vector(strsplit(consensus_data[1,], split = ";")[[1]])
stringr::str_split_fixed(consensus_data$V1, ";", 43) -> data
data <- data[2:nrow(data),]
consensus <- data.frame(data)
colnames(consensus) <- header
data <- consensus

# Extract data, unlisting positions of substitutions
mutation_data_full <- data.frame()
for(r in 1:nrow(data))
{
  
  row <- data[r,]
  id <- unique(row$seqName)
  
  mutations <- str_split(unlist(row[2]), ",")[[1]]
  positions <- as.numeric(substr(mutations, 2, nchar(mutations)-1))
  
  mutation_data <- data.frame(position = positions, id = id, mutation = mutations)
  
  mutation_data_full <- rbind(mutation_data_full, mutation_data)
}

# Combine with metadata
mutation_data_full %>% 
  mutate(ID = gsub("/", "_", id)) %>%
  filter(ID %in% metadata$ID) %>%
  select(-id) -> mutation_data_full

left_join(mutation_data_full, metadata, by = "ID") %>%
  mutate(mutation_RNA = gsub("T", "U", mutation)) -> mutation_data_plot

# Plot

mutation_data_plot_highlightday29 <- filter(mutation_data_plot, mutation_RNA %in% c("C5184U", "C23191U"))
mutation_data_plot_highlightday93 <- filter(mutation_data_plot, mutation_RNA %in% c("C4230U", "C5183U", "A13768C", "C26305U"))
mutation_data_plot_highlightday106 <- filter(mutation_data_plot, mutation_RNA %in% c("C5178U", "C13665U", "C15720U"))

mutation_data_plot$day <- factor(mutation_data_plot$day, levels = rev(c(7, 12, 22, 29, 33, 38, 93, 106)))

genome.substitution.plot <- ggplot() +
  geom_point(data = mutation_data_plot, aes(x = position, y = as.factor(day)), shape = 108, color = "dodgerblue2", size = 7) +
  geom_segment(data = mutation_data_plot, aes(x = 0, y = as.factor(day), xend = 29903, yend = as.factor(day))) +
  ylab("Day of Illness") +
  xlab("Genome Position") +
  theme(axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = c()) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "white")) + theme(legend.position = "")

genome.substitution.plot + 
  geom_point(data = mutation_data_plot_highlightday29, aes(x = position, y = as.factor(day)), color = "lightsalmon1", size = 7, shape = 108) +
  geom_point(data = mutation_data_plot_highlightday93, aes(x = position, y = as.factor(day)), color = "purple", size = 7, shape = 108) +
  geom_point(data = mutation_data_plot_highlightday106, aes(x = position, y = as.factor(day)), color = "firebrick2", size = 7, shape = 108) +
  geom_text(data = mutation_data_plot, aes(x = position, y = as.factor(day), label = mutation_RNA), size = 1, nudge_y = 0.3)

