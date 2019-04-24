Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar");lapply(Packages, library, character.only = TRUE);rm(Packages)
tic()
# Downloaded on 18/03/19 with the options Type=all,organism=human,database=all,evidence=Experimental,confidence=High
RegNet <- read_csv("~/Documents/Thesis/export_Mon_Mar_18_13_56_38_UTC_2019.csv")
dim(RegNet)

# Lets filter the tissue specific interactions
data("rnaGeneTissue")
rnaGeneTissue <- filter(rnaGeneTissue, Sample == "cerebral cortex")
#rnaGeneTissue <- filter(rnaGeneTissue, Sample == "cerebral cortex" & Value != 0)
RegNet_tissue_specific <- semi_join(RegNet, rnaGeneTissue, c("regulator_symbol" = "Gene.name"))
RegNet_tissue_specific <- semi_join(RegNet_tissue_specific, rnaGeneTissue, c("target_symbol" = "Gene.name"))
dim(RegNet_tissue_specific) # 9038    7 and 8439    7 with the Value != 0
write.table(RegNet_tissue_specific[,c(1,3,5)], "/home/tain/Documents/Thesis/RegNet_tissue_specific.txt", sep = "\t", col.names = F, row.names = F)
# Lets run Diamond
#/usr/bin/time -o RegNet_tissue_specific_500.time ./DIAMOND.py RegNet_tissue_specific.txt RBD_list_symbol_only.txt 500 Prediction_RegNet_tissue_specific_500.txt &> RegNet_tissue_specific_500.log
RegNet_tissue_specific_igraph <- igraph::graph_from_data_frame(RegNet_tissue_specific[,c(2,4)], directed = F)
ecount(RegNet_tissue_specific_igraph);vcount(RegNet_tissue_specific_igraph) #5650,2242
RegNet_tissue_specific_EGall <- sort(V(RegNet_tissue_specific_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- RBDlist$X1
names(RBDlist_setA) <- 'setA'

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)
tic()
EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% RegNet_tissue_specific_EGall)
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
RegNet_tissue_specific_neat <- neat::neat(EG2GOBP, RBDlist_setA, RegNet_tissue_specific_igraph, "undirected", RegNet_tissue_specific_EGall, 0.05)
toc()

RegNet_tissue_specific_neat$padj <- p.adjust(RegNet_tissue_specific_neat$pvalue, method = 'BH')
RegNet_tissue_specific_neat_10_most_enriched_GOID <- RegNet_tissue_specific_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)
prediction_RegNet_tissue_specific_500 <- read_delim("~/Documents/Thesis/Prediction_RegNet_tissue_specific_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_RegNet_tissue_specific_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = EG2GOBP[RegNet_tissue_specific_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = RegNet_tissue_specific_igraph, nettype = "undirected", nodes = RegNet_tissue_specific_EGall, alpha = 0.05)
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
