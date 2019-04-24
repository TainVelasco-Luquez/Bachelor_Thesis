tic()
metabosignal_ready <- read_delim("~/Documents/Thesis/metabosignal_ready.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
metabosignal_ready_igraph <- graph_from_data_frame(metabosignal_ready[,1:2])
metabosignal_igraph <- metabosignal_ready_igraph
metabosignal_ready_EGall <- sort(V(metabosignal_ready_igraph)$name)
metabosignal_EGall <- metabosignal_ready_EGall

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
tic()
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% metabosignal_ready_EGall)
toc()

tic()
EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(go_id) %>%
  dplyr::select(gene_id) %>%
  split(f = EG2GOBP$go_id) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(gene_id) %>%
      unlist() %>%
      as.character()})
toc()

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- RBDlist$X1
names(RBDlist_setA) <- 'setA'

tic()
metabosignal_neat <- neat::neat(EG2GOBP, RBDlist_setA, metabosignal_ready_igraph, "undirected", metabosignal_ready_EGall, 0.05)
toc() # 23.453 sec elapsed

#RBDlist$X1[RBDlist$X1 %in% V(metabosignal_ready_igraph)$name]
#plot(induced_subgraph(metabosignal_ready_igraph, as.character(RBDlist$X1[RBDlist$X1 %in% V(metabosignal_ready_igraph)$name])))
# Indeed there are no edges between the nodes
metabosignal_neat$padj <- p.adjust(metabosignal_neat$pvalue, method = 'BH')
metabosignal_neat_10_most_enriched_GOID <- metabosignal_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#/usr/bin/time -o prediction_metabosignal_500.time ./DIAMOND.py metabosignal.txt RBD_list_symbol.txt 500 prediction_metabosignal_500.txt &> prediction_metabosignal_500.log
#/usr/bin/time -o prediction_metabosignal_1000.time ./DIAMOND.py metabosignal.txt RBD_list_symbol.txt 1000 prediction_metabosignal_1000.txt &> prediction_metabosignal_1000.log

prediction_metabosignal_500 <- read_delim("~/Documents/Thesis/prediction_metabosignal_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_metabosignal_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = EG2GOBP[metabosignal_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = metabosignal_igraph, nettype = "undirected", nodes = metabosignal_EGall, alpha = 0.05)
                          }
print(Sys.time() - strt) 
stopCluster(cl) 
validationGOBP$Iteration <- rep(1:500, each = 10)
validationGOBP <- as.data.frame(validationGOBP)
validationGOBP$padj <- p.adjust(validationGOBP$pvalue, method = 'BH')

validationGOBPlot <- ggplot(data = validationGOBP,
                            aes(x = Iteration, y = pvalue, colour = A)) +
  geom_line() +
  scale_y_continuous(trans = "reverse") +
  labs(x = "Iteration", y = "Adjusted NEAT p-value") +
  scale_colour_discrete(name = "GO Biological\nProcess term") +
  geom_hline(aes(yintercept = 0.05),
             color = "tomato",
             linetype = "dashed",
             size = 1,
             show.legend = FALSE) +
  theme(panel.grid.minor = element_blank())
validationGOBPlotly <- ggplotly(validationGOBPlot)
validationGOBPlotly
write.csv(validationGOBP, file = "/home/tain/Documents/Thesis/mwtabosignal_validationGOBP.csv")
toc()

