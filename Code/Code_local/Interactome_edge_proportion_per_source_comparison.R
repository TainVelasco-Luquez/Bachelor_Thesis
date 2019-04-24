#### Packages ####
Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar", "poweRlaw", "ergm", "intergraph", "ggnet", "gridExtra")
lapply(Packages, library, character.only = TRUE)
rm(Packages)

#### Interactome APID 2 + RegNet assembly ####
metabosignal_ready <- read_delim("Data/Individual/metabosignal_ready.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # 27200 edges
names(metabosignal_ready)

APID_level2_ready <- read_delim("Data/Individual/APID_level2_ready.txt", "\t", escape_double = FALSE, col_names = FALSE)
APID_level2_ready$X3 <- "ppi"
names(APID_level2_ready) <- names(metabosignal_ready)

regnet_ready <- read_delim("Data/Individual/regnet_ready.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(regnet_ready) <- names(metabosignal_ready)

interactome_apid2_compounds <- unique(dplyr::bind_rows(mutate_all(APID_level2_ready, as.character), mutate_all(metabosignal_ready, as.character), mutate_all(regnet_ready, as.character)))

metabosignal_ready <- metabosignal_ready %>% filter(!grepl("C", metabosignal_ready$EG_node1), !grepl("C", metabosignal_ready$EG_node2)) # 16116 edges after filtering out those interactions with compounds
interactome_apid2 <- unique(dplyr::bind_rows(
  APID_level2_ready,
  mutate_at(metabosignal_ready, c("EG_node1", "EG_node2"), as.numeric),
  regnet_ready
))
interactome_apid2$Source <- as.factor(interactome_apid2$Source)

#### Comparison of interactomes with differetn sources ####
# Metabosignal (with compounds) + APID level 2 + RegNet
dim(interactome_apid2_compounds) # 51438    3

# Metabosignal (without compounds) + APID level 2 + RegNet
dim(interactome_apid2) # 40354    3

# Metabosignal (without compounds) + APID level 3 + RegNet
interactome_regnet_ready <- read_delim("Data/Individual/interactome_regnet_ready.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# Metabosignal (with compounds) + APID level 3 + Marbach et al (2016)
interactome_ready_2 <- read_delim("Data/Individual/interactome_ready_2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#### Save RDS and importing it ####
save.image(file = "Data/Interactome_edge_proportion_per_source_comparison.RData")
load(file = "Data/Interactome_edge_proportion_per_source_comparison.RData")

#### Function to extract the legend of a ggplot2 figure ####
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### Proportion of interactions ####
interactome_apid2_compounds_barproportion <- ggplot(interactome_apid2_compounds, aes(x = "", fill = Source)) +
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "") +
  ggtitle("Metabosignal (with compounds) + APID level 2 + RegNet")

interactome_apid2_barproportion <- ggplot(interactome_apid2, aes(x = "", fill = Source, stat = "count")) +
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "Proportion of edges") +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  ggtitle("Metabosignal (without compounds) + APID level 2 + RegNet")

interactome_regnet_ready_barproportion <- ggplot(interactome_regnet_ready, aes(x = "", fill = Source)) +
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "Proportion of edges") +
  ggtitle("Metabosignal (without compounds) + APID level 3 + RegNet")

interactome_ready_2_barproportion <- ggplot(interactome_ready_2, aes(x = "", fill = Source)) +
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "") +
  ggtitle("Metabosignal (with compounds) + APID level 3 + Marbach et al (2016)")

my_legend <- g_legend(interactome_apid2_barproportion)
grid.arrange(arrangeGrob(interactome_regnet_ready_barproportion + theme(legend.position = "none"), interactome_apid2_compounds_barproportion + theme(legend.position = "none"), interactome_apid2_barproportion + theme(legend.position = "none"), interactome_ready_2_barproportion + theme(legend.position = "none")), my_legend, heights = c(10, 1))