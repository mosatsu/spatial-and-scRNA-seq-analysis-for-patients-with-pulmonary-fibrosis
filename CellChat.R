# GSE135839 #

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"))
devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
gc()

seurat_object <- readRDS("/GSE135893.rds") # After SCT normalization
data.input <- seurat_object[["SCT"]]$data # normalized data matrix
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels)) 
# Recreate the CellChat object
cellChat <- createCellChat(object = seurat_object, group.by = "celltype", assay = "SCT")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Configure the CellChat database
CellChatDB <- CellChatDB.human
cellChat@DB <- CellChatDB
# Set to NULL if using all features
features <- NULL
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# subset
# cellChat <- subsetData(cellChat, features = features)
future::plan("multisession", workers = 24) # do parallel
options(future.globals.maxSize = 96 * 1024^3)  # 96GB
print(dim(cellChat@data))  
print(rownames(cellChat@data)[1:10])  
# Check if genes and ligand-receptor pairs are properly configured
cellChat <- subsetData(cellChat)  # Subset the ligand-receptor data based on metadata.
# Identify and adjust overexpressed genes
cellChat <- identifyOverExpressedGenes(cellChat)
# Identify overexpressed ligand-receptor pairs
cellChat <- identifyOverExpressedInteractions(cellChat)
# Ensure again that the ligand-receptor pairs are nonzero
if (nrow(cellChat@LR$LRsig) == 0) {
  stop("No highly variable ligand-receptor pairs detected. Please check your data.")
} else {
  print(paste("Number of ligand-receptor pairs: ", nrow(cellChat@LR$LRsig)))
}
ptm = Sys.time()

# Compute cell-cell communication probabilities using the triMean method and stores the result in cellchat.
cellchat <- computeCommunProb(cellChat, type = "triMean")
# Filter the communication data in cellchat, keeping only interactions involving at least 10 cells.
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute communication probabilities at the pathway level and updates cellchat
cellchat <- computeCommunProbPathway(cellchat)
# Aggregate the cell-cell communication network across different signaling pathways in cellchat.
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Step 1: Select the top 10 cell types with the strongest interactions
top10_types <- names(sort(rowSums(cellchat@net$weight), decreasing = TRUE))[1:10]
# Step 2: Restrict the network data to the top 10 cell types
sub_weight <- cellchat@net$weight[top10_types, top10_types]
sub_count <- cellchat@net$count[top10_types, top10_types]
# Filter groupSize to retain only the sizes of relevant cell types
sub_groupSize <- groupSize[match(top10_types, names(groupSize))]
# Step 3: Handle NA values in sub_groupSize (addressing matching failures).
if (any(is.na(sub_groupSize))) {
  warning("NA values found in sub_groupSize! Replacing them with 1.")
  sub_groupSize[is.na(sub_groupSize)] <- 1
}
# Step 4: Visualize the network.
par(mfrow = c(1, 2), xpd = TRUE)

# Visualize the number of interaction
netVisual_circle(sub_count, vertex.weight = sub_groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions (Top 10 Cell Types)")
# Visualize the strength of interaction
netVisual_circle(sub_weight, vertex.weight = sub_groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction weights/strength (Top 10 Cell Types)")

mat <- cellchat@net$weight
par(mfrow = c(1,1), xpd=TRUE)

output_folder <- ".../feature"
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE) 
}
# Save files in a loop
for (i in 1:nrow(mat)) {
  # Replace invalid characters in file names (e.g., /, +)
  file_name <- gsub("/", "_", rownames(mat)[i])  
  file_name <- gsub("\\+", "_plus", file_name)
  file_path <- file.path(output_folder, paste0(file_name, "_circle_plot.pdf"))
  pdf(file = file_path, width = 24, height = 24)
  # Create a temporary matrix and plot it
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])

  dev.off()
}

# Check the available signaling pathways
available_pathways <- cellchat@netP$pathways
print(available_pathways)

# Save
pathways.show <- c("COLLAGEN")
vertex.receiver = seq(1,4) 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
pathways.show <- c("COLLAGEN")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
group.cellType <- c(rep("Myofibroblasts", 14), rep("HAS1 High Fibroblasts", 10), rep("KRT5-/KRT17+", 15)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
pathways.show <- c("COLLAGEN")
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))


# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)
# Conduct CCC analysis for individual ligand-receptor pairs
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


setwd("/Desktop")
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)

#levels(cellchat@idents)
#[1] "AT1"                            "AT2"                            "B Cells"                        "Basal"                          "cDCs"                          
#[6] "Ciliated"                       "Differentiating Ciliated"       "Endothelial Cells"              "Fibroblasts"                    "HAS1 High Fibroblasts"         
#[11] "KRT5-/KRT17+"                   "Lymphatic Endothelial Cells"    "Macrophages"                    "Mast Cells"                     "Mesothelial Cells"             
#[16] "Monocytes"                      "MUC5AC+ High"                   "MUC5B+"                         "Myofibroblasts"                 "NK Cells"                      
#[21] "pDCs"                           "Plasma Cells"                   "PLIN2+ Fibroblasts"             "Proliferating Epithelial Cells" "Proliferating Macrophages"     
#[26] "Proliferating T Cells"          "SCGB3A2+"                       "SCGB3A2+ SCGB1A1+"              "Smooth Muscle Cells"            "T Cells"                       
#[31] "Transitional AT2"   


# Bubble plot
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), remove.isolate = FALSE)
#> Comparing communications on a single object

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), signaling = c("COLLAGEN","TGFb"), remove.isolate = FALSE)
#> Comparing communications on a single object

# (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("COLLAGEN"))
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object


# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
netVisual_bubble(cellchat, targets.use = c("KRT5-/KRT17+"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.target = T)
# (6) sort cell pairs based on the sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T)
# (7) sort cell pairs based on the sources.use and then targets.use defined by users
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T)
# (8) sort cell pairs based on the targets.use and then sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)
# "HAS1 High Fibroblasts"

# Chord diagram
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), lab.cex = 0.5,legend.pos.y = 30)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c("Myofibroblasts"), targets.use = c("KRT5-/KRT17+"), slot.name = "netP", legend.pos.x = 5)
netVisual_chord_gene(cellchat, sources.use = c("HAS1 High Fibroblasts"), targets.use = c("KRT5-/KRT17+"), slot.name = "netP", legend.pos.x = 5)
# Step 1: HAS1 High Fibroblasts と Myofibroblasts → Obtain the signaling intensity of KRT5-/KRT17+
pathway_matrix <- cellchat@netP$prob  # Signaling pathway probability (interaction strength)
# Retrieve the signaling intensity between the specified cell populations
interaction_strengths_by_pathway <- pathway_matrix[c("Myofibroblasts", "HAS1 High Fibroblasts"), 
                                                   "KRT5-/KRT17+", , drop = FALSE]
# Convert to dataframe
interaction_df <- data.frame(
  Pathway = dimnames(interaction_strengths_by_pathway)[[3]],  
  Strength = as.numeric(interaction_strengths_by_pathway)   
)
# Step 2: Rank signaling pathways by interaction strength and select the top 10.
interaction_df <- interaction_df[order(-interaction_df$Strength), ]  
top10_pathways <- interaction_df$Pathway[1:10]  
# Step 3: Display only the selected signaling pathways
netVisual_chord_gene(cellchat, 
                     sources.use = c("Myofibroblasts", "HAS1 High Fibroblasts"), 
                     targets.use = c("KRT5-/KRT17+"), 
                     signaling = top5_pathways,  
                     slot.name = "netP", 
                     legend.pos.x = 10)


plotGeneExpression(cellchat, signaling = "COLLAGEN", enriched.only = TRUE, type = "violin")

# USERS can show the expression of all signaling genes related to one signaling pathway by
plotGeneExpression(cellchat, signaling = "COLLAGEN", enriched.only = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

## additional analysis ##
# Obtain the signaling intensities of all signaling pathways from HAS1 High Fibroblasts to KRT5-/KRT17+
pathway_matrix <- cellchat@netP$prob  # Interaction strength of each signaling pathway
# Retrieve interaction strength for each signaling pathway
interaction_strengths_by_pathway <- pathway_matrix[source_cell, target_cell, , drop = FALSE]
str(interaction_strengths_by_pathway)
# Convert to dataframe
interaction_df <- data.frame(
  Pathway = dimnames(interaction_strengths_by_pathway)[[3]],  
  Strength = as.numeric(interaction_strengths_by_pathway)     
)
interaction_df <- interaction_df[order(-interaction_df$Strength), ]

print(interaction_df)

