
# Pre-filtered Visium-seurat data made from the h5 file and spatial data
fHP <- readRDS("fHP.rds")
CTDILD <- readRDS("CTDILD.rds")
IPF <- readRDS("IPF.rds")
uILD <- readRDS("uILD.rds")

fHP$orig.ident <- "fHP"
CTDILD$orig.ident <- "CTDILD"
IPF$orig.ident <- "IPF"
uILD$orig.ident <- "uILD"

merged <- merge(
  x = fHP,
  y = list(CTDILD, IPF, uILD),
  add.cell.ids = c("fHP", "CTDILD", "IPF", "uILD")
)

# ---- SCTransform normalization ---- #
options(future.globals.maxSize = 96 * 1024^3)  # 96GB
merged <- SCTransform(merged, assay = "Spatial", method = "glmGamPoi", verbose = FALSE)


SM2_genes <- c("FKBP10", "DCN", "COL5A2", "COL6A1", "CERCAM", "FBN1",
               "TNC", "COL16A1", "THBS2", "MMP2", "COL5A1", "POSTN",
               "CTHRC1", "COL6A2", "LUM", "COL6A3", "SPARC", "COL1A1",
               "COL3A1", "COL1A2")

merged <- AddModuleScore(
  object = merged,
  features = list(SM2_genes),
  name = "SM2_Score",
  assay = "SCT"  
)
SpatialFeaturePlot(merged, features = "SM2_Score1", min.cutoff = "0") 


# ---- Integration ---- #
features <- SelectIntegrationFeatures(object.list = list(fHP, CTDILD, IPF, uILD), nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = list(fHP, CTDILD, IPF, uILD), anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)
sample_info <- merged$orig.ident
names(sample_info) <- colnames(integrated)
integrated <- AddMetaData(integrated, metadata = sample_info, col.name = "orig.ident")
DefaultAssay(integrated) <- "integrated"  # "integrated"アッセイを使用
integrated <- ScaleData(integrated, features = VariableFeatures(integrated), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
# ---- Harmony ---- #
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident")
# ---- UMAP---- #
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)
# A list of genes present in the current object
all_genes <- rownames(integrated)
# Check the presence of SM2_genes
valid_genes <- SM2_genes[SM2_genes %in% all_genes]
# Check the absence of SM2_genes
missing_genes <- SM2_genes[!SM2_genes %in% all_genes]
print(valid_genes)   
print(missing_genes) 
integrated <- AddModuleScore(
  object = integrated,
  features = list(SM2_genes),
  name = "SM2_Score",
  assay = "integrated"
)
SpatialFeaturePlot(integrated, features = "SM2_Score1", min.cutoff = "0", max.cutoff = "40") 


# ---- Calculate the proportion of high-score regions ---- #
# median or preconfigured thrashold
threshold <- median(integrated$SM2_Score1)
threshold <- 20
  # Convert orig.ident to a factor type and specify its order
      integrated@meta.data$orig.ident <- factor(
        integrated@meta.data$orig.ident,
        levels = c("fHP", "CTDILD", "IPF", "uILD")  
      )
      
      # Calculate the proportion of high-score spots for each ILD pattern
      high_score_stats <- integrated@meta.data %>%
        group_by(orig.ident) %>%
        summarise(
          Total_Spots = n(),
          High_Score_Spots = sum(SM2_Score1 >= threshold),
          High_Score_Percentage = (High_Score_Spots / Total_Spots) * 100
        )
      # Results
      print(high_score_stats)
      # ----Statistical analysis among ILD patterns ---- #
      kruskal_test <- kruskal.test(SM2_Score1 ~ orig.ident, data = integrated@meta.data)
      print(kruskal_test)
      # ---- Barplot ---- #
      ggplot(high_score_stats, aes(x = orig.ident, y = High_Score_Percentage, fill = orig.ident)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        ylab("Percentage of High SM2 Score Spots") +
        xlab("ILD Pattern") +
        ggtitle("Comparison of High SM2 Score Regions across ILD Patterns")

      ggplot(high_score_stats, aes(x = orig.ident, y = High_Score_Percentage, fill = orig.ident)) +
        geom_bar(stat = "identity") +
        theme_minimal(base_size = 14) +  
        ylab("Percentage of High SM2 Score Spots") +
        xlab("ILD Pattern") +
        ggtitle("Comparison of High SM2 Score Regions across ILD Patterns") +
        theme(
          axis.text = element_text(size = 14),       
          axis.title = element_text(size = 14),      
          plot.title = element_text(size = 14, hjust = 0.5),  
          legend.title = element_text(size = 14),    
          legend.text = element_text(size = 12)   
        )
      