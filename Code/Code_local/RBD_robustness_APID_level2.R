#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
Packages <- c("ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly", "hpar"); lapply(Packages, library, character.only = TRUE); rm(Packages)
tic() # 39.699 sec elapsed
set.seed(123)

APID_level2 <- read_delim("~/Documents/Thesis/APID_level2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
APID_level2 <- APID_level2[, c(1,2,5)] #171448 edges
APID_level2_igraph <- graph_from_data_frame(APID_level2[,2:3])
vcount(APID_level2_igraph)
ecount(APID_level2_igraph)
is.simple(APID_level2_igraph)
APID_level2_EGall <- sort(V(APID_level2_igraph)$name)
UNIPROTID2ENSEMBL <- AnnotationDbi::select(org.Hs.eg.db, APID_level2_EGall, c("ENSEMBL", "UNIPROT"), "UNIPROT")

APID_level2_ENSEMBLID <- dplyr::inner_join(APID_level2, UNIPROTID2ENSEMBL, c("UniprotID_A" = "UNIPROT")); APID_level2_ENSEMBLID <- dplyr::inner_join(APID_level2_ENSEMBLID, UNIPROTID2ENSEMBL, c("UniprotID_B" = "UNIPROT"))
# Because UNIPROTID to ENSEMBLID has 1 to many mapping, number of edges in the network increase. This is an artifact as it is the same vertex but with different ID, so it is required to remove this. Fortunatelly, APID has an unique interactionID:
APID_level2_ENSEMBLID <- dplyr::distinct(APID_level2_ENSEMBLID, InteractionID, .keep_all = T) 
dim(APID_level2_ENSEMBLID) # 171448 edges. However some of them have NAs.  Lets drop such edges:
any(is.na(APID_level2_ENSEMBLID))
APID_level2_ENSEMBLID <- APID_level2_ENSEMBLID[complete.cases(APID_level2_ENSEMBLID),]
dim(APID_level2_ENSEMBLID) # 135069 edges left after the ID mapping

data("hpaNormalTissue")
hpaNormalTissue <- dplyr::filter(hpaNormalTissue,
                                 Reliability == "Supportive" | Reliability == "Approved",
                                 Level == "High" | Level == "Low" | Level == "Medium",
                                 Tissue == "hippocampus" | Tissue == "hypothalamus" | Tissue == "caudate" | Tissue == "cerebellum" | Tissue == "cerebral cortex")

APID_level2_ready <- dplyr::semi_join(APID_level2_ENSEMBLID,
                                      hpaNormalTissue,
                                      by = c("ENSEMBL.x" = "Gene"));APID_level2_ready <- dplyr::semi_join(APID_level2_ready,
                                      hpaNormalTissue,
                                      by = c("ENSEMBL.y" = "Gene"))
dim(APID_level2_ready) # 15364 edges left after tissue filtering
APID_level2_ready_igraph <- graph_from_data_frame(APID_level2_ready[,4:5])
ecount(APID_level2_ready_igraph) # 15364
vcount(APID_level2_ready_igraph) # 3740

APID_level2_ready_EGall <- sort(V(APID_level2_ready_igraph)$name)
ENSEMBL2EG <- AnnotationDbi::select(org.Hs.eg.db, APID_level2_ready_EGall, c("ENSEMBL", "ENTREZID"), "ENSEMBL")
APID_level2_ready_EGID <- dplyr::inner_join(APID_level2_ready, ENSEMBL2EG, c("ENSEMBL.x" = "ENSEMBL"))
APID_level2_ready_EGID <- dplyr::inner_join(APID_level2_ready_EGID, ENSEMBL2EG, c("ENSEMBL.y" = "ENSEMBL"))
APID_level2_ready_EGID <- dplyr::distinct(APID_level2_ready_EGID, InteractionID, .keep_all = T) 

APID_level2_ready_EGID <- APID_level2_ready_EGID[,6:7]
APID_level2_ready_EGID <- APID_level2_ready_EGID[complete.cases(APID_level2_ready_EGID),]
APID_level2_ready_EGID$ENTREZID.x <- as.factor(APID_level2_ready_EGID$ENTREZID.x)
APID_level2_ready_EGID$ENTREZID.y <- as.factor(APID_level2_ready_EGID$ENTREZID.y)
write.table(APID_level2_ready_EGID, file = "/home/tain/Documents/Thesis/APID_level2_ready.txt", sep = "\t", col.names = F, quote = F, row.names = F)
APID_level2_ready_igraph <- igraph::graph_from_data_frame(APID_level2_ready_EGID, directed = F)
ecount(APID_level2_ready_igraph);vcount(APID_level2_ready_igraph) #15364,3739
APID_level2_ready_EGall <- sort(V(APID_level2_ready_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- RBDlist$X1
names(RBDlist_setA) <- 'setA'

org.Hs.egGODB <- AnnotationDbi::as.data.frame(org.Hs.egGO)

EG2GOBP <- dplyr::filter(org.Hs.egGODB, Ontology == "BP", Evidence == "EXP" | Evidence == "IDA" | Evidence == "IMP" | Evidence == "IGI" | Evidence == "IEP" | Evidence == "ISS" | Evidence == "ISA" | Evidence == "ISO") %>% dplyr::arrange(go_id) %>% dplyr::select(go_id, gene_id) %>% unique() %>% dplyr::filter(gene_id %in% APID_level2_ready_EGall)

EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(go_id) %>%
  dplyr::select(gene_id) %>%
  split(f = EG2GOBP$go_id) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(gene_id) %>%
      unlist() %>%
      as.character()})

tic()
APID_level2_ready_neat <- neat::neat(EG2GOBP, RBDlist_setA, APID_level2_ready_igraph, "undirected", APID_level2_ready_EGall, 0.05)
toc() # 9.065 sec elapsed

RBDlist$X1[RBDlist$X1 %in% V(APID_level2_ready_igraph)$name]
plot(induced_subgraph(APID_level2_ready_igraph, as.character(RBDlist$X1[RBDlist$X1 %in% V(APID_level2_ready_igraph)$name])))
# Indeed there are no edges between the nodes

APID_level2_ready_neat$padj <- p.adjust(APID_level2_ready_neat$pvalue, method = 'BH')
APID_level2_ready_neat_10_most_enriched_GOID <- APID_level2_ready_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#/usr/bin/time -o prediction_APID_level2_ready_500.time ./DIAMOND.py APID_level2_ready.txt RBD_list_symbol.txt 500 prediction_APID_level2_ready_500.txt &> prediction_APID_level2_ready_500.log
#/usr/bin/time -o prediction_APID_level2_ready_1000.time ./DIAMOND.py APID_level2_ready.txt RBD_list_symbol.txt 1000 prediction_APID_level2_ready_1000.txt &> prediction_APID_level2_ready_1000.log

prediction_APID_level2_ready_500 <- read_delim("~/Documents/Thesis/prediction_APID_level2_ready_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_APID_level2_ready_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = EG2GOBP[APID_level2_ready_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = APID_level2_ready_igraph, nettype = "undirected", nodes = APID_level2_ready_EGall, alpha = 0.05)
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
write.csv(validationGOBP, file = "/home/tain/Documents/Thesis/APID_level2_validationGOBP.csv")
toc()