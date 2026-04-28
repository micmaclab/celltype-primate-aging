library(Matrix)
library(Seurat)
library(CellChat)
library(patchwork)


run_cellchat <- function(cellchat) {
  cellchat@idents <- droplevels(cellchat@idents)
  
  # 1. Data Preparation
  message(paste0("[", Sys.time(), "] Step 1: Subsetting and identifying overexpressed genes..."))
  cellchat <- subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 2. Probability Computation 
  message(paste0("[", Sys.time(), "] Step 3: Computing communication probabilities (this may take a while)..."))
  cellchat <- computeCommunProb(
    cellchat, raw.use = TRUE, 
    type ='truncatedMean', trim = 0.05
  )
  
  # 3. Filtering
  message(paste0("[", Sys.time(), "] Step 4: Filtering communication (min.cells = 10)..."))
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # 4. Pathway and Network Aggregation
  message(paste0("[", Sys.time(), "] Step 5: Computing pathway probabilities and aggregating networks..."))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  message(paste0("[", Sys.time(), "] Done! CellChat object is ready for visualization."))
  return(cellchat)
}
CellChatDB <- CellChatDB.human 

# Create V1 cellchat
spa_matrix_V1 <- readRDS("spa_matrix_V1.rds")

#Remove missing values
spa_matrix_V1 <- subset(
  spa_matrix_V1,
  subset = !is.na(SubClass) & SubClass != "NA" & SubClass != ""
)

# Set correct level of cell class as the identity
spa_matrix_V1$SubClass <- droplevels(as.factor(spa_matrix_V1$SubClass))
Idents(spa_matrix_V1) <- "SubClass"

# downsample each subclass up to 50,000 cells
spa_matrix_V1 <- subset(spa_matrix_V1, downsample = 50000)
print(ncol(spa_matrix_V1))

# create cell chat object
cellchat_V1 <- createCellChat(
  object = spa_matrix_V1,
  group.by = "SubClass"
)

# set database of ligand/receptor pairs to compare with
cellchat_V1@DB <- CellChatDB

# run cell chat function above
cellchat_V1 <-run_cellchat(cellchat_V1)

# save cellchat 
file_name <- "cellchat_V1_results_raw_data_trim0.05.rds"
saveRDS(cellchat_V1, file_name)
message(paste0("Saved: ",file_name))


# load in dlpfc seurat object
spa_matrix_dlpfc <- readRDS("spa_matrix_dlpfc.rds")

# remove missing labels
spa_matrix_dlpfc <- subset(
  spa_matrix_dlpfc,
  subset = !is.na(SubClass) & SubClass != "NA" & SubClass != ""
)

# reset cell identities to align with subclass 
spa_matrix_dlpfc$SubClass <- droplevels(as.factor(spa_matrix_dlpfc$SubClass))
Idents(spa_matrix_dlpfc) <- "SubClass"

# downsample 500,000 cells per subclass
spa_matrix_dlpfc <- subset(spa_matrix_dlpfc, downsample = 50000)
print(ncol(spa_matrix_dlpfc))

# create cellchat object
cellchat_dlpfc <- createCellChat(
  object = spa_matrix_dlpfc,
  group.by = "SubClass"
)

# set database of ligand/receptor pairs to compare with
cellchat_dlpfc@DB <- CellChatDB

# run cellchat code above
cellchat_dlpfc <-run_cellchat(cellchat_dlpfc)

# save cellchat object
file_name <- "cellchat_dlpfc_results_raw_data_trim0.05.rds"
saveRDS(cellchat_dlpfc, file_name)
message(paste0("Saved: ",file_name))

