Packages <- c("ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar"); lapply(Packages, library, character.only = TRUE); rm(Packages)
brain_adult <- read_delim("~/Documents/Thesis/brain_adult.csv",
                          "\t", escape_double = FALSE, col_names = c("Node1", "Node2", "Weight"),
                          trim_ws = TRUE)
brain_adult_igraph <- graph_from_data_frame(brain_adult[,1:2])
# Lets check if there is a natural cutoff for cl or LCC in the brain adult network to choose as the threshold to include in the network
source("/home/tain/MyRfunctions/my_targeted_robustness.R")
tic();robustness_brain_adult_igraph <- my_targeted_robustness(brain_adult_igraph, 22);toc()
# lets plot the results
plot(robustness_brain_adult_igraph$PercentRemoved, robustness_brain_adult_igraph$cl)
plot(robustness_brain_adult_igraph$PercentRemoved, robustness_brain_adult_igraph$LCC)
# There is a cutoff when aproximately 10 percent of data is removed, which is still to high as this will include too many interaftions compared againste metabosignal and PPII
# Lets do the same procedure for the scaling parameter 
load("~/Documents/Thesis/Bachelor_Thesis/Code/Data/HIvalidation.RData")
odin <- degree(interactomeGraph)
freya <- fit_power_law(odin, implementation = "plfit") # alpha=3.7, xmin=2920
tic();loki <- my_alpha_on_deletion(brain_adult_igraph, 22);toc()
loki <- loki[complete.cases(loki),]
plot(loki$PercentRetained, loki$alpha)
