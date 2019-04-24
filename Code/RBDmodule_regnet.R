#### Packages ####
Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar", "poweRlaw", "ergm", "intergraph", "ggnet", "GO.db"); lapply(Packages, library, character.only = TRUE); rm(Packages)

#### Reading in files ####
interactome <- read_delim("~/Documents/Thesis/interactome_regnet_ready.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
interactome$Source <- as.factor(interactome$Source)
dim(interactome) # 30074     3

#### Proportion of interactions ####
int1barproportion <- ggplot(interactome, aes(x = "", fill = Source)) +
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "Proportion of edges") +
  ggtitle("A")
int1barproportion

#### Graph features ####
interactomeGraph <- graph_from_data_frame(interactome, directed = F)
summary(interactomeGraph)
# Saving the number of supported sources as a weight attribute
E(interactomeGraph)$Weight <- igraph::count.multiple(interactomeGraph)
# Removing loops
interactomeGraph <- igraph::simplify(interactomeGraph, remove.loops = TRUE, remove.multiple = FALSE)
# Extracting topological properties of interest
vcount(interactomeGraph)
ecount(interactomeGraph)
mean(degree(interactomeGraph, normalized = FALSE))
median(degree(interactomeGraph, normalized = FALSE))
max(degree(interactomeGraph, normalized = FALSE))
min(degree(interactomeGraph, normalized = FALSE))
transitivity(interactomeGraph)
average.path.length(interactomeGraph)
diameter(interactomeGraph)

# Checking number of components and selecting only the LCC
table(sapply(decompose.graph(interactomeGraph), vcount))
interactomeGraph <- decompose.graph(interactomeGraph)[[1]]

#### Plotting the distribution ####
CumDegDis <- ggplot() +
  geom_point(aes(
    c(1:length(degree.distribution(interactomeGraph,
      cumulative = TRUE
    ))),
    degree.distribution(interactomeGraph,
      cumulative = TRUE
    )
  ),
  colour = "black",
  size = 1
  ) +
  geom_vline(aes(
    xintercept = mean(igraph::degree(interactomeGraph,
      normalized = FALSE
    )),
    color = "Mean"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = T
  ) +
  geom_vline(aes(
    xintercept = median(igraph::degree(interactomeGraph,
      normalized = FALSE
    )),
    color = "Median"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = T
  ) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Degree (k)", y = "P(k)") +
  scale_color_manual(name = "Statistics", values = c("Mean" = "red", "Median" = "green")) +
  ggtitle("A") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
CumDegDis
####  Distribution per source ####
# $P(k)_{metabosignal}$
metabosignalGraph <- graph_from_data_frame(interactome[interactome$Source == "metabosignal", ])
CumDegDisMetaboSignal <- ggplot() +
  geom_point(aes(
    c(1:length(degree.distribution(metabosignalGraph,
      cumulative = TRUE
    ))),
    degree.distribution(metabosignalGraph,
      cumulative = TRUE
    )
  ),
  colour = "tomato",
  size = 1
  ) +
  geom_vline(aes(
    xintercept = mean(igraph::degree(metabosignalGraph,
      normalized = FALSE
    )),
    color = "Mean"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  geom_vline(aes(
    xintercept = median(igraph::degree(metabosignalGraph,
      normalized = FALSE
    )),
    color = "Median"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Degree (k)", y = "P(k)") +
  scale_color_manual(name = "Statistics", values = c("Mean" = "red", "Median" = "green")) +
  ggtitle("C") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
CumDegDisMetaboSignal

# $P(k)_{PPI}$
ppiGraph <- graph_from_data_frame(interactome[interactome$Source == "ppi", ])
CumDegDisPPI <- ggplot() +
  geom_point(aes(
    c(1:length(degree.distribution(ppiGraph,
      cumulative = TRUE
    ))),
    degree.distribution(ppiGraph,
      cumulative = TRUE
    )
  ),
  colour = "gold",
  size = 1
  ) +
  geom_vline(aes(
    xintercept = mean(igraph::degree(ppiGraph,
      normalized = FALSE
    )),
    color = "Mean"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  geom_vline(aes(
    xintercept = median(igraph::degree(ppiGraph,
      normalized = FALSE
    )),
    color = "Median"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Degree (k)", y = "P(k)") +
  scale_color_manual(name = "Statistics", values = c("Mean" = "red", "Median" = "green")) +
  ggtitle("D") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
CumDegDisPPI

# $P(k)_{regulatory}$
regulatoryGraph <- graph_from_data_frame(interactome[interactome$Source == "regulatory", ])
CumDegDisRegulatory <- ggplot() +
  geom_point(aes(
    c(1:length(degree.distribution(regulatoryGraph,
      cumulative = TRUE
    ))),
    degree.distribution(regulatoryGraph,
      cumulative = TRUE
    )
  ),
  colour = "gray50",
  size = 1
  ) +
  geom_vline(aes(
    xintercept = mean(igraph::degree(regulatoryGraph,
      normalized = FALSE
    )),
    color = "Mean"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  geom_vline(aes(
    xintercept = median(igraph::degree(regulatoryGraph,
      normalized = FALSE
    )),
    color = "Median"
  ),
  linetype = "dashed",
  size = 1,
  show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Degree (k)", y = "P(k)") +
  scale_color_manual(name = "Statistics", values = c("Mean" = "red", "Median" = "green")) +
  ggtitle("B") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
CumDegDisRegulatory
#### Scale free? ####
# Extracting the degree of interactomeGraph
degInteractome <- igraph::degree(interactomeGraph,
  normalized = FALSE,
  v = V(interactomeGraph)
)
interactomeGraphPower <- displ$new(degInteractome)
interactomeGraphPower$setXmin(estimate_xmin(interactomeGraphPower,
  xmins = NULL,
  pars = NULL,
  xmax = max(igraph::degree(interactomeGraph,
    normalized = FALSE
  )),
  distance = "ks"
))
interactomeGraphPower$getPars()
interactomeGraphPower$getXmin()
get_ntail(interactomeGraphPower, prop = FALSE, lower = FALSE)

tic()
bootsinteractomeGraphPower <- bootstrap_p(interactomeGraphPower, no_of_sims = 1000, threads = detectCores() - 1, seed = 0)
toc() # 105.477 sec elapsed
bootsinteractomeGraphPower$p # 0.921
bootsinteractomeGraphPower$gof # 0.04891

#### KNN ####
interactomeGraphKNN <- interactomeGraph
E(interactomeGraphKNN)$weight <- 1
interactomeGraphKNNWeighted <- simplify(interactomeGraphKNN)
E(interactomeGraphKNNWeighted)$weight
interactomeGraphKNN <- knn(interactomeGraphKNNWeighted, vids = V(interactomeGraphKNNWeighted), weights = NULL)
plot(interactomeGraphKNN$knn, log = "xy", col = "goldenrod", xlab = c("Log-Degree"), ylab = c("Log Average Neighbor Degree"))

average.path.length(interactomeGraph)
diameter(interactomeGraph)
transitivity(interactomeGraph)

## Assortativity
assortativity_degree(RBDmoduleGraph)

#### HI validation ####
library(MyRfunctions)
tic()
interactome_cl <- MyRfunctions::cl_against_random(interactomeGraph, 1000, 25)
toc() # 6824.293 sec elapsed

HItrans.plot <- ggplot(data = interactome_cl$simulated_df, aes(x = Transitivity)) +
  geom_histogram(binwidth = 0.0008, alpha = .5) +
  geom_vline(aes(xintercept = mean(interactome_cl$simulated_df[1:1000, 1], na.rm = TRUE)), linetype = "dashed", size = 0.7, color = "red") +
  ylab("P(cl)") +
  xlab("cl") +
  geom_segment(aes(x = interactome_cl$observed_cl, y = 200, xend = interactome_cl$observed_cl, yend = 10),
    arrow = arrow(length = unit(0.3, "cm")), colour = "#00BFC4", alpha = 0.8, size = 0.8
  ) +
  annotate("text", x = interactome_cl$observed_cl, y = 280, label = paste("Observed\ncl = ", round(interactome_cl$observed_cl, digits = 3)), size = 4)
HItrans.plot

tic()
interactome_l <- MyRfunctions::l_against_random(interactomeGraph, 1000, 23)
toc() # 146.293 sec elapsed

HIapl.plot <- ggplot(data = interactome_l$simulated_df, aes(x = APL)) +
  geom_histogram(binwidth = 0.0008, alpha = .5) +
  geom_vline(aes(xintercept = mean(interactome_l$simulated_df[1:1000, 1], na.rm = TRUE)), linetype = "dashed", size = 0.7, color = "red") +
  ylab("P(< l >)") +
  xlab("l") +
  geom_segment(aes(x = interactome_l$observed_l, y = 18, xend = interactome_l$observed_l, yend = 2),
    arrow = arrow(length = unit(0.3, "cm")), colour = "#00BFC4", alpha = 0.8, size = 0.8
  ) +
  annotate("text", x = interactome_l$observed_l, y = 25, label = paste("Observed\nl = ", round(interactome_l$observed_l, digits = 3)), size = 4)
HIapl.plot

#### HI ERGM fitting ####
interactomeGraph.net <- intergraph::asNetwork(interactomeGraph)
interactomeGraph.net
tic()
interactomeGraph.model <- ergm(interactomeGraph.net ~ edges, control = control.ergm(parallel = 23, parallel.type = "PSOCK"))
toc() # 7.057 sec elapsed
summary(interactomeGraph.model)
# Log odds of any tie derived from ergm?
plogis(coef(interactomeGraph.model)[["edges"]]) # 0.0018846 it is very unlikely that HI ties are random
tic()
interactomeGraph.model.gof.deg <- gof(interactomeGraph.model ~ degree, control = control.gof.ergm(nsim = 100, parallel = 25, parallel.type = "PSOCK"), verbose = T)
toc() # 2438.545 sec elapsed and 2700.565 sec elapsed
plot(interactomeGraph.model.gof.deg)
plot(interactomeGraph.model.gof.deg, plotlogodds = TRUE, xlim = c(0, 30))

#### RBD Module ####
RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist <- as.character(RBDlist$X1)
# filtering out nodes not present in the interactome
RBDlist_in_HI <- V(interactomeGraph)$name[V(interactomeGraph)$name %in% RBDlist]

#### RBD Proto-Module ####
protoModule <- igraph::induced_subgraph(graph = interactomeGraph, vids = RBDlist_in_HI, impl = "auto")

# Converting from igraph to netork class, for visualisation with ggnet2
protoModuleNet <- intergraph::asNetwork(protoModule)

# Plotting disconnected nodes
protoModulePlot <- ggnet2(protoModuleNet,
  mode = "circle",
  color = "#00BFC4",
  alpha = 0.5,
  edge.alpha = 0.9,
  edge.size = 0.1,
  label = TRUE,
  label.size = 5,
  label.alpha = 1
) +
  geom_point(aes(color = color), size = 5, alpha = 0.6) +
  ggtitle("A")
protoModulePlot

#### DIAMOnD ####

# Lets run Diamond
# /usr/bin/time -o interactome_RegNet_500.time ./DIAMOND.py interactome_regnet_ready.txt RBD_list.txt 500 Prediction_interactome_RegNet_500.txt &> interactome_RegNet_500.log

#### Seed nodes NEAT ####

interactome_igraph <- igraph::graph_from_data_frame(interactome, directed = F)
ecount(interactome_igraph)
vcount(interactome_igraph) # 5650,2242
interactome_EGall <- sort(V(interactome_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector("list", 1)
RBDlist_setA[[1]] <- RBDlist$X1
names(RBDlist_setA) <- "setA"

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% interactome_EGall)

EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(go_id) %>%
  dplyr::select(gene_id) %>%
  split(f = EG2GOBP$go_id) %>%
  lapply(FUN = function(i) {
    as.data.frame(i) %>%
      dplyr::select(gene_id) %>%
      unlist() %>%
      as.character()
  })

tic()
interactome_neat <- neat::neat(EG2GOBP, RBDlist_setA, interactome_igraph, "undirected", interactome_EGall, 0.05)
toc() # 20.975 sec elapsed

interactome_neat$padj <- p.adjust(interactome_neat$pvalue, method = "BH")
interactome_neat_10_most_enriched_GOID <- interactome_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#### DIAMOnD module NEAT ####

prediction_interactome_500 <- read_delim("~/Documents/Thesis/Prediction_interactome_RegNet_500.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)
setB <- vector("list", 1)
setB[[1]] <- prediction_interactome_500$DIAMOnD_node
names(setB) <- "setB"

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl)
strt <- Sys.time()
validationGOBP <- foreach(
  i = 1:length(setB$setB),
  .combine = "rbind"
) %dopar% {
  RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
  neat::neat(alist = EG2GOBP[interactome_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = interactome_igraph, nettype = "undirected", nodes = interactome_EGall, alpha = 0.05)
}
print(Sys.time() - strt) # 20.096 mins
stopCluster(cl)
validationGOBP$Iteration <- rep(1:500, each = 10)
validationGOBP <- as.data.frame(validationGOBP)
validationGOBP$padj <- p.adjust(validationGOBP$pvalue, method = "BH")

validationGOBPlot <- ggplot(
  data = validationGOBP,
  aes(x = Iteration, y = padj, colour = A)
) +
  geom_line() +
  scale_y_continuous(trans = "reverse") +
  labs(x = "Iteration", y = "Adjusted NEAT p-value") +
  scale_colour_discrete(name = "GO Biological\nProcess term") +
  geom_hline(aes(yintercept = 0.05),
    color = "tomato",
    linetype = "dashed",
    size = 1,
    show.legend = FALSE
  ) +
  theme(panel.grid.minor = element_blank())
validationGOBPlotly <- ggplotly(validationGOBPlot)
validationGOBPlotly

write.table(validationGOBP, file = "~/Documents/Thesis/RegNet.csv", sep = ",", col.names = T, quote = F, row.names = F)

#### Module plot ####
Prediction_interactome_RegNet_500 <- read_delim("~/Documents/Thesis/Prediction_interactome_RegNet_500.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(Prediction_interactome_RegNet_500) <- c("rank", "name")
Prediction_interactome_RegNet_500 <- Prediction_interactome_RegNet_500[1:405,2] 
Prediction_interactome_RegNet_500$Source <- "DIAMOnD"

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist$Source <- "Seed"
names(RBDlist) <- c("name", "Source")

RBDmodule <- bind_rows(RBDlist, Prediction_interactome_RegNet_500)
RBDmodule$name <- as.character(unlist(RBDmodule$name))
ENTREZ2SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(RBDmodule$name)), c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME"), "ENTREZID")
RBDmodule <- dplyr::left_join(RBDmodule, ENTREZ2SYMBOL, c("name" = "ENTREZID"))
RBDmodule <- dplyr::distinct(RBDmodule, name, .keep_all = T) 

write.table(RBDmodule, file = "RBD_module.txt", quote = F, row.names = F, col.names = T, sep = "\t")

vids <- V(interactomeGraph)$name[V(interactomeGraph)$name %in% RBDmodule$name]

RBDmoduleGraph <- igraph::induced_subgraph(graph = interactomeGraph, vids = vids, impl = "auto")
RBDmoduleGraph <- igraph::simplify(RBDmoduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

RBDmoduleNet <- asNetwork(RBDmoduleGraph)

RBDmoduleNet %v% "Degree" <- sna::degree(RBDmoduleNet)
RBDmoduleNet %v% "Type" <- ifelse(get.vertex.attribute(RBDmoduleNet, "vertex.names") %in% as.character(unlist(RBDlist$name)), "Seed", "DIAMOnD")

#RBDmoduleNet %e% "Type" <- ifelse(get.vertex.attribute(RBDmoduleNet, "vertex.names") %in% RBDlist$Source == "Seed", "#00BFC4", "#F8766D")

RBD_list_symbol <- read_delim("~/Documents/Thesis/RBD_list_symbol.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# Plotting disconnected nodes
RBDmoduleNetPlot_layout <- gplot.layout.fruchtermanreingold(RBDmoduleNet, NULL)
RBDmoduleNet %v% "x" <- RBDmoduleNetPlot_layout[, 1]
RBDmoduleNet %v% "y" <- RBDmoduleNetPlot_layout[, 2]

RBDmoduleNetPlot <- ggnet2(RBDmoduleNet,
  mode = c("x", "y"),
  label = odin$X3,
  node.size = "Degree",
  node.color = "Type",
  alpha = 0.5,
  edge.alpha = 0.7,
  edge.size = 0.1,
  palette = c("Seed" = "#00BFC4", "DIAMOnD" = "#F8766D"),
  label.size = 3,
  legend.position = "none"
) +
  geom_point(aes(color = color), color = "white", alpha = 0.2, show.legend = FALSE) +
  geom_point(aes(color = color), alpha = 0.5)
RBDmoduleNetPlot

write.table(get.data.frame(RBDmoduleGraph), file = "RBDmoduleGraph.txt", sep = "\t", quote = F, row.names = F)

#### Module validation ####
tic()
RBDmodule_cl <- MyRfunctions::cl_against_random(RBDmoduleGraph, 1000, 25)
toc() # 2.315 sec elapsed

RBDmodule_trans.plot <- ggplot(data = RBDmodule_cl$simulated_df, aes(x = Transitivity)) +
  geom_histogram(binwidth = 0.0008, alpha = .5) +
  geom_vline(aes(xintercept = mean(RBDmodule_cl$simulated_df[1:1000, 1], na.rm = TRUE)), linetype = "dashed", size = 0.7, color = "red") +
  ylab("P(cl)") +
  xlab("cl") +
  geom_segment(aes(x = RBDmodule_cl$observed_cl, y = 200, xend = RBDmodule_cl$observed_cl, yend = 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "#00BFC4", alpha = 0.8, size = 0.8
  ) +
  annotate("text", x = RBDmodule_cl$observed_cl, y = 280, label = paste("Observed\ncl = ", round(RBDmodule_cl$observed_cl, digits = 3)), size = 4);RBDmodule_trans.plot


tic()
RBDmodule_l <- MyRfunctions::l_against_random(RBDmoduleGraph, 1000, 23)
toc() # 2.536 sec elapsed

RBDmodule_apl.plot <- ggplot(data = RBDmodule_l$simulated_df, aes(x = APL)) +
  geom_histogram(binwidth = 0.0008, alpha = .5) +
  geom_vline(aes(xintercept = mean(RBDmodule_l$simulated_df[1:1000, 1], na.rm = TRUE)), linetype = "dashed", size = 0.7, color = "red") +
  ylab("P(< l >)") +
  xlab("l") +
  geom_segment(aes(x = RBDmodule_l$observed_l, y = 18, xend = RBDmodule_l$observed_l, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "#00BFC4", alpha = 0.8, size = 0.8
  ) +
  annotate("text", x = RBDmodule_l$observed_l, y = 25, label = paste("Observed\nl = ", round(RBDmodule_l$observed_l, digits = 3)), size = 4);RBDmodule_apl.plot

#### Module ERGM fitting ####
RBDmoduleGraph.net <- intergraph::asNetwork(RBDmoduleGraph)
RBDmoduleGraph.net
tic()
RBDmoduleGraph.model <- ergm(RBDmoduleGraph.net ~ edges, control = control.ergm(parallel = 23, parallel.type = "PSOCK"))
toc() # 1.085 sec elapsed
summary(RBDmoduleGraph.model)
# Log odds of any tie derived from ergm?
plogis(coef(RBDmoduleGraph.model)[["edges"]]) # 0.031237 it is unlikely that HI ties are random
tic()
RBDmoduleGraph.model.gof.deg <- gof(RBDmoduleGraph.model ~ degree, control = control.gof.ergm(nsim = 100, parallel = 25, parallel.type = "PSOCK"), verbose = T)
toc() # 1603.933 sec elapsed
par(mfrow = c(1,2))
plot(RBDmoduleGraph.model.gof.deg)
plot(RBDmoduleGraph.model.gof.deg, plotlogodds = TRUE, xlim = c(0, 30))
par(mfrow = c(1,1))

#### Module modularity ####
RBDmoduleGraph_layout <- layout_with_fr(RBDmoduleGraph, niter = 10000)

RBDmodule_partition <- fastgreedy.community(RBDmoduleGraph)
sizes(RBDmodule_partition)
membership(RBDmodule_partition)
modularity(RBDmodule_partition)
plot(RBDmodule_partition, RBDmoduleGraph, layout = RBDmoduleGraph_layout, vertex.label = NA, vertex.size = 5)

RBDmodule_laplacian <- graph.laplacian(RBDmoduleGraph)
RBDmodule_eigen <- eigen(RBDmodule_laplacian)
plot(RBDmodule_eigen$values, col="blue", ylab="Eigenvalues of Graph Laplacian")

RBDmodule_walktrap <- walktrap.community(RBDmoduleGraph)
sizes(RBDmodule_walktrap)
plot(RBDmodule_walktrap, RBDmoduleGraph, vertex.label = NA)

RBDmodule_betweeness.community <- edge.betweenness.community(RBDmoduleGraph)
sizes(RBDmodule_betweeness.community)

#### Module GOBP enrichment ####
RBDmoduleGraph <- read_delim("RBDmoduleGraph.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
RBDmoduleGraph <- graph_from_data_frame(RBDmoduleGraph, directed = F)
RBDmodule_EGall <- sort(V(RBDmoduleGraph)$name)

RBDmodule <- read_delim("~/Documents/Thesis/RBDmodule.csv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
RBDmodule_setA <- vector("list", 1)
RBDmodule_setA[[1]] <- RBDmodule$X1
names(RBDmodule_setA) <- "setA"

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% RBDmodule_EGall)
EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(go_id) %>%
  dplyr::select(gene_id) %>%
  split(f = EG2GOBP$go_id) %>%
  lapply(FUN = function(i) {
    as.data.frame(i) %>%  
      dplyr::select(gene_id) %>%
      unlist() %>%
      as.character()
  })

# Remeber to set the network and the nodes IDs to the interactome instead of the Module.
tic()
RBDmodulerich <- neat::neat(EG2GOBP, RBDmodule_setA, interactomeGraph, "undirected", interactome_EGall, 0.05)
toc() #1.227 sec elapsed

RBDmodulerich$padj <- p.adjust(RBDmodulerich$p, method = 'BH')
RBDmodulerich <- arrange(as_tibble(RBDmodulerich),  desc(padj))
summary(RBDmodulerich)
GOID2GOTERM <- AnnotationDbi::select(GO.db, as.character(unique(unlist(RBDmodulerich$A))), c("DEFINITION", "TERM", "ONTOLOGY"), "GOID")
RBDmodulerich <- left_join(RBDmodulerich, GOID2GOTERM, c("A" = "GOID"))
  
write.table(RBDmodulerich, file = "~/Documents/Thesis/RBDmodulerich.csv", sep = "\t", col.names = T, quote = F, row.names = F)

#Top 10 most significantly enriched GOIDs 
RBDmodulerich10 <- RBDmodulerich %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)
 
#### Error handling #### 
# A very common error of NETA is: "There are no edges connected to genes in GO:0001921 . The test cannot be computed." Which appears when the query set is not a vector inside a list but instead another class inside the list. e.g. (see the "RBD Module" section for context) RBDlist <- as.character(RBDlist) assigns the whole dataframe that only contains 1 vector, but it is a DF!. The correct line is RBDlist <- as.character(RBDlist$X1), assigning explicitly only the vector.
