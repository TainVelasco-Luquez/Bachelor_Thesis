Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly")
lapply(Packages, library, character.only = TRUE)
tic() # 39.699 sec elapsed
set.seed(123)
brain_adult <- read_delim("~/Documents/Thesis/brain_adult.csv",
                          "\t", escape_double = FALSE, col_names = c("Node1", "Node2", "Weight"),
                          trim_ws = TRUE)
nrow(brain_adult) * 0.025 # Taking only the 2.5% of interactions sums up 29032 edges
brain_adult_2.5 <- brain_adult %>% arrange(desc(Weight)) %>% top_n(29032)
#write.table(brain_adult_2.5[,1:2], file = "brain_adult_2.5.txt", sep = "\t", col.names = F)
brain_adult_2.5_igraph <- igraph::graph_from_data_frame(brain_adult_2.5, directed = F)
brain_adult_2.5_EGall <- sort(V(brain_adult_2.5_igraph)$name)
SYMBOL2EG <- AnnotationDbi::select(org.Hs.eg.db, brain_adult_2.5_EGall, c("SYMBOL", "ENTREZID"), "SYMBOL")
tic();brain_adult_2.5_EGID <- dplyr::left_join(brain_adult_2.5, SYMBOL2EG, c("Node1" = "SYMBOL"));toc() # 0.011 sec elapsed
tic();brain_adult_2.5_EGID <- dplyr::left_join(brain_adult_2.5_EGID, SYMBOL2EG, c("Node2" = "SYMBOL"));toc() # 0.013 sec elapsed
brain_adult_2.5_EGID <- brain_adult_2.5_EGID[,4:5]
brain_adult_2.5_EGID <- brain_adult_2.5_EGID[complete.cases(brain_adult_2.5_EGID),]
brain_adult_2.5_EGID$ENTREZID.x <- as.factor(brain_adult_2.5_EGID$ENTREZID.x)
brain_adult_2.5_EGID$ENTREZID.y <- as.factor(brain_adult_2.5_EGID$ENTREZID.y)
brain_adult_2.5_igraph <- igraph::graph_from_data_frame(brain_adult_2.5_EGID, directed = F)
brain_adult_2.5_EGall <- sort(V(brain_adult_2.5_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
EZ2SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(RBDlist)), c("SYMBOL", "ENTREZID"), "ENTREZID")
#write.table(EZ2SYMBOL, file="RBD_list_symbol.txt", sep= "\t", col.names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- EZ2SYMBOL$ENTREZID
names(RBDlist_setA) <- 'setA'

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
tic()
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% brain_adult_2.5_EGall)
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

tic()
brain_adult_2.5_neat <- neat::neat(EG2GOBP, RBDlist_setA, brain_adult_2.5_igraph, "undirected", brain_adult_2.5_EGall, 0.05)
toc() # 23.453 sec elapsed
brain_adult_2.5_neat$padj <- p.adjust(brain_adult_2.5_neat$pvalue, method = 'BH')
brain_adult_2.5_neat_10_most_enriched_GOID <- brain_adult_2.5_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#/usr/bin/time -o prediction_brain_adult_2.5_500.time ./DIAMOND.py brain_adult_2.5.txt RBD_list_symbol.txt 500 prediction_brain_adult_2.5_500.txt &> prediction_brain_adult_2.5_500.log
#/usr/bin/time -o prediction_brain_adult_2.5_1000.time ./DIAMOND.py brain_adult_2.5.txt RBD_list_symbol.txt 1000 prediction_brain_adult_2.5_1000.txt &> prediction_brain_adult_2.5_1000.log

prediction_brain_adult_2.5_500 <- read_delim("~/Documents/Thesis/prediction_brain_adult_2.5_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_brain_adult_2.5_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = EG2GOBP[brain_adult_2.5_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = brain_adult_2.5_igraph, nettype = "undirected", nodes = brain_adult_2.5_EGall, alpha = 0.05)
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
toc()