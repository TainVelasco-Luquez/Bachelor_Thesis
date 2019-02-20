Packages <- c("ggplot2", "dplyr", "ggplot2", "org.Hs.eg.db", "readr", "neat", "doParallel", "igraph", "tictoc", "plotly")
lapply(Packages, library, character.only = TRUE)

######### Interactome
interactome_ready <- read_delim("~/Documents/Thesis/interactome_ready_2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(interactome_ready)
dim(interactome_ready[interactome_ready$Source == "metabosignal",]) #27171
dim(interactome_ready[interactome_ready$Source == "regulatory",]) #3024355 
dim(interactome_ready[interactome_ready$Source == "ppi",]) #5020 

# Reducing the number of regulatory interactions by filtering out those whose edge weight is lower than 
brain_adult <- read_delim("~/Documents/Thesis/brain_adult.csv", "\t", escape_double = FALSE, col_names = c("Node1", "Node2", "Weight"),trim_ws = TRUE)
summary(brain_adult$Weight)
sum(brain_adult$Weight >= 0.01725) # Number of edges whose wheight is greater than the mean = 268355
nrow(brain_adult) * 0.10 # Taking only the 10% of interactions sums up 116129 edges
nrow(brain_adult) * 0.05 # Taking only the 5% of interactions sums up 58064 edges
nrow(brain_adult) * 0.025 # Taking only the 2.5% of interactions sums up 29032 edges which is almost the same as the metabosignal inetractions

######## brain_adult 2.5
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
### 1000 iterations
prediction_brain_adult_2.5_3000 <- read_delim("~/Documents/Thesis/prediction_brain_adult_2.5_3000.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_brain_adult_2.5_3000[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = SYMBOL2GOBP[brain_adult_2.5_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = igraph::graph_from_data_frame(brain_adult_2.5), nettype = "undirected", nodes = brain_adult_2.5_EGall, alpha = 0.05)
                          }
print(Sys.time() - strt) 
stopCluster(cl) 
validationGOBP$Iteration <- rep(1:3000, each = 10)
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

######## brain_adult 5
brain_adult <- read_delim("~/Documents/Thesis/brain_adult.csv",
                          "\t", escape_double = FALSE, col_names = c("Node1", "Node2", "Weight"),
                          trim_ws = TRUE)
nrow(brain_adult) * 0.5 # Taking only the 5% of interactions sums up 580643 edges
brain_adult_5 <- brain_adult %>% arrange(desc(Weight)) %>% top_n(580643)
write.table(brain_adult_5[,1:2], file = "brain_adult_5.txt", sep = "\t", col.names = F)
brain_adult_5_igraph <- igraph::graph_from_data_frame(brain_adult_5, directed = F)
brain_adult_5_EGall <- sort(V(brain_adult_5_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
EZ2SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(RBDlist)), c("SYMBOL", "ENTREZID"), "ENTREZID")
write.table(EZ2SYMBOL, file="RBD_list_symbol.txt", sep= "\t", col.names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- EZ2SYMBOL$SYMBOL
names(RBDlist_setA) <- 'setA'

SYMBOL2GOBP <- AnnotationDbi::select(org.Hs.eg.db, brain_adult_5_EGall, c("GO", "SYMBOL"), "SYMBOL")

tic()
SYMBOL2GOBP <- SYMBOL2GOBP %>%
  dplyr::group_by(GO) %>%
  # dplyr::select(ENTREZID) %>%
  split(f = SYMBOL2GOBP$GO) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(SYMBOL) %>%
      unlist() %>%
      as.character()})
toc()

tic()
brain_adult_5_neat <- neat::neat(SYMBOL2GOBP, RBDlist_setA, graph_from_data_frame(brain_adult_5), "undirected", brain_adult_5_EGall, 0.05)
toc() # 23.453 sec elapsed
brain_adult_5_neat$padj <- p.adjust(brain_adult_5_neat$pvalue, method = 'BH')
brain_adult_5_neat_10_most_enriched_GOID <- brain_adult_5_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#/usr/bin/time -o prediction_brain_adult_5_500.time ./DIAMOND.py brain_adult_5.txt RBD_list_symbol.txt 500 prediction_brain_adult_5_500.txt &> prediction_brain_adult_5_500.log

prediction_brain_adult_5_500 <- read_delim("~/Documents/Thesis/prediction_brain_adult_5_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_brain_adult_5_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = SYMBOL2GOBP[brain_adult_5_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = igraph::graph_from_data_frame(brain_adult_5), nettype = "undirected", nodes = brain_adult_5_EGall, alpha = 0.05)
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

######## brain_adult greater than mean
brain_adult <- read_delim("~/Documents/Thesis/brain_adult.csv",
                          "\t", escape_double = FALSE, col_names = c("Node1", "Node2", "Weight"),
                          trim_ws = TRUE)
nrow(brain_adult[brain_adult >= 0.01725]) # Taking only those taht are greater than the mean summing up 268355 edges
brain_adult_5 <- brain_adult %>% arrange(desc(Weight)) %>% top_n(580643)
write.table(brain_adult_5[,1:2], file = "brain_adult_5.txt", sep = "\t", col.names = F)
brain_adult_5_igraph <- igraph::graph_from_data_frame(brain_adult_5, directed = F)
brain_adult_5_EGall <- sort(V(brain_adult_5_igraph)$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
EZ2SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(RBDlist)), c("SYMBOL", "ENTREZID"), "ENTREZID")
write.table(EZ2SYMBOL, file="RBD_list_symbol.txt", sep= "\t", col.names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- EZ2SYMBOL$SYMBOL
names(RBDlist_setA) <- 'setA'

SYMBOL2GOBP <- AnnotationDbi::select(org.Hs.eg.db, brain_adult_5_EGall, c("GO", "SYMBOL"), "SYMBOL")

tic()
SYMBOL2GOBP <- SYMBOL2GOBP %>%
  dplyr::group_by(GO) %>%
  # dplyr::select(ENTREZID) %>%
  split(f = SYMBOL2GOBP$GO) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(SYMBOL) %>%
      unlist() %>%
      as.character()})
toc()

tic()
brain_adult_5_neat <- neat::neat(SYMBOL2GOBP, RBDlist_setA, graph_from_data_frame(brain_adult_5), "undirected", brain_adult_5_EGall, 0.05)
toc() # 23.453 sec elapsed
brain_adult_5_neat$padj <- p.adjust(brain_adult_5_neat$pvalue, method = 'BH')
brain_adult_5_neat_10_most_enriched_GOID <- brain_adult_5_neat %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

#/usr/bin/time -o prediction_brain_adult_5_500.time ./DIAMOND.py brain_adult_5.txt RBD_list_symbol.txt 500 prediction_brain_adult_5_500.txt &> prediction_brain_adult_5_500.log

prediction_brain_adult_5_500 <- read_delim("~/Documents/Thesis/prediction_brain_adult_5_500.txt", "\t",escape_double = FALSE,trim_ws = TRUE, col_names = T)
setB <- vector('list', 1)
setB[[1]] <- unlist(prediction_brain_adult_5_500[, 2], use.names = FALSE)
names(setB) <- 'setB'

cl <- makeCluster(detectCores() - 5)
registerDoParallel(cl) 
strt <- Sys.time() 
validationGOBP <- foreach(i = 1:length(setB$setB),
                          .combine = 'rbind') %dopar% {
                            RBDlist_setA[[1]] <- c(RBDlist_setA[[1]], i)
                            neat::neat(alist = SYMBOL2GOBP[brain_adult_5_neat_10_most_enriched_GOID[, 1] ], blist = RBDlist_setA, network = igraph::graph_from_data_frame(brain_adult_5), nettype = "undirected", nodes = brain_adult_5_EGall, alpha = 0.05)
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

######## Regulatory
regulatory_ready_2 <- read_delim("~/Documents/Thesis/regulatory_ready_2.txt", "\t", escape_double = FALSE)
regulatory_ready_2_EGall <- sort(V(graph_from_data_frame(regulatory_ready_2))$name)
RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
EZ2SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(RBDlist)), c("SYMBOL", "ENTREZID"), "ENTREZID")
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- EZ2SYMBOL$SYMBOL
names(RBDlist_setA) <- 'setA'
EG2GOBP <- AnnotationDbi::select(org.Hs.eg.db, regulatory_ready_2_EGall, c("GO", "ENTREZID"), "ENTREZID")
tic()
EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(GO) %>%
  # dplyr::select(ENTREZID) %>%
  split(f = EG2GOBP$GO) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(ENTREZID) %>%
      unlist() %>%
      as.character()})
toc() # 36.504 sec elapsed
tic()
RBDrich <- neat(alist = EG2GOBP, blist = RBDlist_setA, network = graph_from_data_frame(regulatory_ready_2), nettype = "undirected", nodes = regulatory_ready_2_EGall, alpha = 0.01)
toc()

RBDrich$padj <- p.adjust(RBDrich$pvalue, method = 'BH')
RBDrich10 <- RBDrich %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)
cl <- makeCluster(detectCores() - 3)
registerDoParallel(cl)  # Required for foreach
strt <- Sys.time()  # For controlling the running time
validationGOBP <- foreach(i = 1:length(prediction_ppi_500_setB$prediction_ppi_500_setB),
                          .combine = 'rbind') %dopar% {
                            setA[[1]] <- c(setA[[1]], i)
                            neat::neat(alist = EG2GOBP[RBDrich10[, 1] ], blist = setA, network = graph_from_data_frame(ppi_ready), nettype = "undirected", nodes = prediction_ppi_500_EGall, alpha = 0.05)
                          }
print(Sys.time() - strt)  # Total run time
stopCluster(cl)  # Giving back the resources to the OS
validationGOBP$Iteration <- rep(1:500, each = 10)

####### PPI only
# /usr/bin/time -o prediction_metabosignal_500.time ./DIAMOND.py metabosignal_ready.txt RBD_list.txt 500 prediction_metabosignal_500.txt &> prediction_metabosignal_500.log
# /usr/bin/time -o prediction_ppi_500.time ./DIAMOND.py ppi_ready.txt RBD_list.txt 500 prediction_ppi_500.txt &> prediction_ppi_500.log

ppi_ready <- read_delim("~/Documents/Thesis/ppi_ready.txt", "\t", escape_double = FALSE)
ppi_ready_EGall <- sort(V(graph_from_data_frame(ppi_ready))$name)

RBDlist <- read_delim("~/Documents/Thesis/RBD_list.txt", "\t", escape_double = FALSE, col_names = F)
RBDlist_setA <- vector('list', 1)
RBDlist_setA[[1]] <- RBDlist
names(RBDlist_setA) <- 'setA'

prediction_ppi_500 <- read_delim("~/Documents/Thesis/prediction_ppi_500.txt",
"\t", escape_double = FALSE, col_names = TRUE,trim_ws = TRUE)
prediction_ppi_500_setB <- vector('list', 1)
prediction_ppi_500_setB[[1]] <- unlist(prediction_ppi_500[, 2], use.names = FALSE)
names(prediction_ppi_500_setB) <- 'setB'

library("org.Hs.eg.db", lib.loc="~/Rmylib")
EG2GOBP <- select(org.Hs.eg.db, ppi_ready_EGall, c("GO", "ENTREZID"), "ENTREZID")

EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(GO) %>%
  dplyr::select(ENTREZID) %>%
  split(f = EG2GOBP$GO) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(ENTREZID) %>%
      unlist() %>%
      as.character()})

RBDrich <- neat(alist = EG2GOBP, blist = RBDlist_setA, network = graph_from_data_frame(ppi_ready), nettype = "undirected", nodes = ppi_ready_EGall, alpha = 0.01)
RBDrich$padj <- p.adjust(RBDrich$pvalue, method = 'BH')
RBDrich10 <- RBDrich %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)

library("doParallel", lib.loc="~/Rmylib")
cl <- makeCluster(detectCores() - 3)
registerDoParallel(cl)  # Required for foreach
strt <- Sys.time()  # For controlling the running time
validationGOBP <- foreach(i = 1:length(prediction_ppi_500_setB$prediction_ppi_500_setB),
                          .combine = 'rbind') %dopar% {
                            setA[[1]] <- c(setA[[1]], i)
                            neat::neat(alist = EG2GOBP[RBDrich10[, 1] ], blist = setA, network = graph_from_data_frame(ppi_ready), nettype = "undirected", nodes = prediction_ppi_500_EGall, alpha = 0.05)
                          }
print(Sys.time() - strt)  # Total run time
stopCluster(cl)  # Giving back the resources to the OS
validationGOBP$Iteration <- rep(1:500, each = 10)

######## Metabosignal
metabosignal_ready <- read_delim("~/Documents/Thesis/metabosignal_ready.txt", "\t", escape_double = FALSE)
metabosignal_ready_EGall <- sort(V(graph_from_data_frame(metabosignal_ready))$name)

prediction_metabosignal_500 <- read_delim("~/Documents/Thesis/prediction_metabosignal_500.txt",
                                 "\t", escape_double = FALSE, col_names = TRUE,trim_ws = TRUE)
prediction_metabosignal_500_setB <- vector('list', 1)
prediction_metabosignal_500_setB[[1]] <- unlist(prediction_metabosignal_500[, 2], use.names = FALSE)
names(prediction_metabosignal_500_setB) <- 'setB'

library("org.Hs.eg.db", lib.loc="~/Rmylib")
EG2GOBP <- select(org.Hs.eg.db, metabosignal_ready_EGall, c("GO", "ENTREZID"), "ENTREZID")

EG2GOBP <- EG2GOBP %>%
  dplyr::group_by(GO) %>%
  dplyr::select(ENTREZID) %>%
  split(f = EG2GOBP$GO) %>%
  lapply(FUN = function(i){ 
    as.data.frame(i) %>% 
      dplyr::select(ENTREZID) %>%
      unlist() %>%
      as.character()})

metabosignal_RBDrich <- neat(alist = EG2GOBP, blist = RBDlist_setA, network = graph_from_data_frame(metabosignal_ready), nettype = "undirected", nodes = metabosignal_ready_EGall, alpha = 0.01)
RBDrich$padj <- p.adjust(RBDrich$paslue, method = 'BH')

prediction_metabosignal_500 <- read_delim("~/Documents/Thesis/prediction_metabosignal_500.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
prediction_metabosignal_500_setB <- vector('list', 1)
prediction_metabosignal_500_setB[[1]] <- unlist(prediction_metabosignal_500[, 2], use.names = FALSE)
names(prediction_metabosignal_500_setB) <- 'setB'

#####
SYMBOL2EG <- AnnotationDbi::select(org.Hs.eg.db, brain_adult_2.5_EGall, c("SYMBOL", "ENTREZID"), "SYMBOL")
tic();brain_adult_2.5_EGID <- dplyr::left_join(brain_adult_2.5, SYMBOL2EG, c("Node1" = "SYMBOL"));toc() # 0.011 sec elapsed
tic();brain_adult_2.5_EGID <- dplyr::left_join(brain_adult_2.5_EGID, SYMBOL2EG, c("Node2" = "SYMBOL"));toc() # 0.013 sec elapsed
length(unique(brain_adult_2.5$Node1)) # Number of features 584
length(unique(brain_adult_2.5$Node2))  # Number of features 3631
length(unique(brain_adult_2.5_EGID$ENTREZID.x)) # After ID conversion are 543 left (41 lost)
length(unique(brain_adult_2.5_EGID$ENTREZID.y)) # After ID conversion are 3477 left (154 lost)
anyNA(brain_adult_2.5)
anyNA(brain_adult_2.5_EGID)
sum(is.na(brain_adult_2.5_EGID$ENTREZID.y))
sum(is.na(brain_adult_2.5_EGID$ENTREZID.x))

brain_adult_2.5_EGID <- brain_adult_2.5_EGID[,4:5]
brain_adult_2.5_EGID <- brain_adult_2.5_EGID[complete.cases(brain_adult_2.5_EGID),]

length(unique(brain_adult_2.5_EGID$ENTREZID.x)) # After remove NA are 541 left (2 lost)
length(unique(brain_adult_2.5_EGID$ENTREZID.y)) # After remove NA are 3443left (34 lost)

sum(is.na(brain_adult_2.5_EGID$ENTREZID.x))
sum(is.na(brain_adult_2.5_EGID$ENTREZID.y))

brain_adult_2.5_EGID$ENTREZID.x <- as.factor(brain_adult_2.5_EGID$ENTREZID.x)
brain_adult_2.5_EGID$ENTREZID.y <- as.factor(brain_adult_2.5_EGID$ENTREZID.y)


#######
alist <- EG2GOBP
for (i in 1:length(alist)) {
  if (is.factor(alist[[i]]) == T) {
    alist[[i]] = as.character(alist[[i]])
  }
  eina[[i]] = (net[, 1] %in% alist[[i]]) | (net[, 
                                                2] %in% alist[[i]])
  if (sum(eina[[i]]) == 1) 
    netred[[i]] = as.matrix(t(net[eina[[i]], ]))
  else if (sum(eina[[i]]) == 0) 
    stop(paste("There are no edges connected to genes in", 
               names(alist)[i], ". The test cannot be computed."))
  else netred[[i]] = net[eina[[i]], ]
  alogic[[i]] = (netred[[i]][, 1] %in% alist[[i]])
  alogic2[[i]] = (netred[[i]][, 2] %in% alist[[i]])
  oa[i] = sum(alogic[[i]]) + sum(alogic2[[i]])
}

odin <- lapply(EG2GOBP, function(i){ifelse(RBDlist_setA$setA %in% i, TRUE, FALSE)})

#########

