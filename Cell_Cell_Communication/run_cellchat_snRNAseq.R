

library(Matrix)
library(Seurat)
library(CellChat)
library(patchwork)

run_cellchat <- function(cellchat) {
  
  # 1. Data Preparation
  message(paste0("[", Sys.time(), "] Step 1: Subsetting and identifying overexpressed genes..."))
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 2. Probability Computation 
  message(paste0("[", Sys.time(), "] Step 2: Computing communication probabilities (this may take a while)..."))
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean") 
  
  # 3. Filtering
  message(paste0("[", Sys.time(), "] Step 3: Filtering communication (min.cells = 10)..."))
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # 4. Pathway and Network Aggregation
  message(paste0("[", Sys.time(), "] Step 4: Computing pathway probabilities and aggregating networks..."))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  message(paste0("[", Sys.time(), "] Done! CellChat object is ready for visualization."))
  return(cellchat)
}

spa_matrix_m1 <- readRDS("spa_matrix_m1.rds")
cellchat_m1 <-run_cellchat(spa_matrix_m1)
saveRDS(cellchat_m1, "cellchat_m1.rds")

spa_matrix_m2 <- readRDS("spa_matrix_m2.rds")
cellchat_m2 <-run_cellchat(spa_matrix_m2)
saveRDS(cellchat_m2, "cellchat_m2.rds")