Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar");lapply(Packages, library, character.only = TRUE);rm(Packages)
tic()
# Downloaded on 18/03/19 with the options Type=all,organism=human,database=all,evidence=Experimental,confidence=High
interactome <- read_delim("~/Documents/Thesis/interactome_regnet_ready.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
dim(interactome) # 9726    7

# Lets run Diamond
#/usr/bin/time -o interactome_RegNet_500.time ./DIAMOND.py interactome_regnet_ready.txt RBD_list.txt 500 Prediction_interactome_RegNet_500.txt &> interactome_RegNet_500.log

interactome_igraph <- igraph::graph_from_data_frame(interactome, directed = F)
ecount(interactome_igraph);vcount(interactome_igraph) #5650,2242
interactome_EGall <- sort(V(interactome_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- RBDlist$X1
names(RBDlist_setA) <- 'setA'

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
tic()
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% interactome_EGall)
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
interactome_neat <- neat::neat(EG2GOBP, RBDlist_setA, interactome_igraph, "undirected", interactome_EGall, 0.05)
toc()

interactome_neat$padj <- p.adjust(interactome_neat$pvalue, method = 'BH')
interactome_neat_10_most_enriched_GOID <- interactome_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)
prediction_interactome_500 <- read_delim("~/Documents/Thesis/Prediction_interactome_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_interactome_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = EG2GOBP[interactome_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = interactome_igraph, nettype = "undirected", nodes = interactome_EGall, alpha = 0.05)
                          }
print(Sys.time() - strt) 
stopCluster(cl) 
validationGOBP$Iteration <- rep(1:500, each = 10)
validationGOBP <- as.data.frame(validationGOBP)
validationGOBP$padj <- p.adjust(validationGOBP$pvalue, method = 'BH')

validationGOBPlot <- ggplot(data = validationGOBP,
                            aes(x = Iteration, y = padj, colour = A)) +
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

write.table(validationGOBP, file = "~/Documents/Thesis/RegNet.csv", sep = ",", col.names = T, quote = F, row.names = F)
toc()
