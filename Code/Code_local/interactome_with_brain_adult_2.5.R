brain_adult_2_5 <- read_delim("~/Documents/Thesis/brain_adult_2.5.txt",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(brain_adult_2_5) <- c("InteractionID", "Node1", "Node2")

keys <- unique(brain_adult_2_5$Node1); keys <- c(keys, unique(brain_adult_2_5$Node2))
SYMBOL2EG <- AnnotationDbi::select(org.Hs.eg.db, keys, c("ENTREZID", "SYMBOL"), "SYMBOL")
brain_adult_2_5 <- dplyr::left_join(brain_adult_2_5, SYMBOL2EG, c("Node1" = "SYMBOL")); brain_adult_2_5 <- dplyr::left_join(brain_adult_2_5, SYMBOL2EG, c("Node2" = "SYMBOL"))
brain_adult_2_5 <- brain_adult_2_5[,4:5]
brain_adult_2_5$Source <- "regulatory"
names(brain_adult_2_5) <- c("EG_node1", "EG_node2", "Source")

metabosignal_ready <- read_delim("~/Documents/Thesis/metabosignal_ready.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)

ppi_ready <- read_delim("~/Documents/Thesis/ppi_ready.txt", "\t", escape_double = FALSE)
ppi_ready$EG_node1 <- as.character(ppi_ready$EG_node1); ppi_ready$EG_node2 <- as.character(ppi_ready$EG_node2)

interactome <- unique(dplyr::bind_rows(ppi_ready,
                                       metabosignal_ready,
                                       brain_adult_2_5))

ggplot(interactome ,aes(x = "", fill = Source)) + 
  geom_bar(position = "fill", width = 0.3) +
  scale_fill_manual(values = c("tomato", "gold", "gray50")) +
  labs(x = "", y = "Proportion of edges") +
  ggtitle("A") 

interactomeGraph <- graph_from_data_frame(interactome)
mean(degree(interactomeGraph, normalized = FALSE))
median(degree(interactomeGraph, normalized = FALSE))
max(degree(interactomeGraph, normalized = FALSE))
min(degree(interactomeGraph, normalized = FALSE))

table(sapply(decompose.graph(interactomeGraph), vcount))

ggplot() +
  geom_point(aes(c(1:length(degree.distribution(interactomeGraph,
                                                cumulative = TRUE)
  )),
  degree.distribution(interactomeGraph,
                      cumulative = TRUE)),
  colour = 'black',
  size = 1) +
  geom_vline(aes(xintercept = mean(igraph::degree(interactomeGraph,
                                                  normalized = FALSE)),
                 color = "Mean"),
             linetype = "dashed",
             size = 1,
             show.legend = T) +
  geom_vline(aes(xintercept = median(igraph::degree(interactomeGraph,
                                                    normalized = FALSE)),
                 color = "Median"),
             linetype = "dashed",
             size = 1,
             show.legend = T) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Degree (k)", y = "P(k)") +
  scale_color_manual(name = "Statistics", values = c("Mean" = "red", "Median" = "green")) +
  ggtitle("A") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
