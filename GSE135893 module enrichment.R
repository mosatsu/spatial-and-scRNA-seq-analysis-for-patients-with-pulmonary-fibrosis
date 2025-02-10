## Make Seurat object and conduct the enrichment analysis ##
gc()
# Define the file path
mtx_path <- ".../GSE135893_matrix.mtx"
barcode_path <- ".../GSE135893_barcodes.tsv"
gene_path <- ".../GSE135893_genes.tsv"
# load
matrix <- Matrix::readMM(mtx_path)
barcodes <- readLines(barcode_path)
genes <- readLines(gene_path)
# Assign barcodes and gene names to the matrix
colnames(matrix) <- barcodes
rownames(matrix) <- genes
head(matrix) 

# load sample meta data #
meta_data_path <- "..../GSE135893_IPF_metadata.csv"
meta_data <- read.table(meta_data_path, header=TRUE, sep=",", stringsAsFactors=FALSE, row.names=1)

# integrate matrix into meta-data
# The count of cell barcode are difference between in matrix and in meta-data
common_cells <- intersect(colnames(matrix), rownames(meta_data))
matrix_filtered <- matrix[, common_cells]
meta_data_filtered <- meta_data[common_cells, ]
# Make SingleCellExperiment object
sce <- SingleCellExperiment(assays=list(counts=matrix_filtered), colData=meta_data_filtered)

# get counts data and meta-data from sce object
sce_counts <- counts(sce)
sce_meta_data <- as.data.frame(colData(sce))

# Make Seurat object
seurat_obj <- CreateSeuratObject(counts = sce_counts, 
                                 min.cells = 3, 
                                 min.features = 200, 
                                 meta.data = sce_meta_data)

#### QC ####
# Add mitochondrial RNA sequences that begin with "MT-" to the data as the "percent.mt" column
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# low-quality cell filtering
seurat_obj <- subset(seurat_obj, subset = 
                       percent.mt < 15 & 
                       nFeature_RNA < 8000)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#### Normalization・scaling ####
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
options(future.globals.maxSize = 96 * 1024^3)  # 96GBに設定（必要に応じて増減）
seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(seurat_obj, ".../GSE135893.rds")
gc()
seurat_obj <- readRDS(".../GSE135893.rds")
seurat_obj <- subset(seurat_obj, subset = Diagnosis %in% c("cHP", "Control", "ILD", "IPF", "NSIP"))
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
DimPlot(seurat_obj, reduction = "umap", group.by = "Diagnosis", pt.size = 0.1, label.size = 5, label=FALSE, raster=FALSE)
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", pt.size = 0.1, label.size = 5, label=FALSE, raster=FALSE)
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", pt.size = 0.1, label.size = 1, label=TRUE, raster=FALSE)
seurat_obj <- AddModuleScore(object = seurat_obj,
                       features = list(c("SPARC", "COL1A1", "MMP2", "THBS2", "CERCAM", "FBN1", "TNC", "POSTN", "LUM", "FKBP10", "CTHRC1", 
                                         "COL3A1", "COL1A2", "COL6A3", "COL6A2", "COL5A1", "COL16A1", "COL5A2", "DCN","COL6A1")),
                       name = "SM2_score1")
FeaturePlot(seurat_obj,"SM2_score11", min.cutoff = 0)
ILDs <- subset(seurat_obj, subset = celltype %in% c("Myofibroblasts", "Fibroblasts", "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts"))
ILDs <- RunPCA(ILDs)
ILDs <- FindNeighbors(ILDs, dims = 1:20)
ILDs <- FindClusters(ILDs, resolution = 0.5)
ILDs <- RunUMAP(ILDs, dims = 1:20)
DimPlot(ILDs, reduction = "umap", group.by = "Diagnosis", pt.size = 0.5, label.size = 5, label=FALSE, raster=FALSE)
DimPlot(ILDs, reduction = "umap", group.by = "celltype", pt.size = 0.5, label.size = 5, label=FALSE, raster=FALSE)
ILDs <- AddModuleScore(object = ILDs,
                             features = list(c("SPARC", "COL1A1", "MMP2", "THBS2", "CERCAM", "FBN1", "TNC", "POSTN", "LUM", "FKBP10", "CTHRC1", 
                                               "COL3A1", "COL1A2", "COL6A3", "COL6A2", "COL5A1", "COL16A1", "COL5A2", "DCN","COL6A1")),
                             name = "SM2_score1")
FeaturePlot(ILDs,"SM2_score11", min.cutoff = 0.75)
FeaturePlot(ILDs,"SPARC", min.cutoff = 2)
FeaturePlot(ILDs,"COL1A1", min.cutoff = 2)
FeaturePlot(ILDs,"COL1A2", min.cutoff = 2)
FeaturePlot(ILDs,"COL3A1", min.cutoff = 2)
DotPlot(ILDs, features = c("COL1A1","COL1A2","COL3A1", "SPARC"), group.by = "celltype", scale = TRUE) +coord_flip()
DotPlot(ILDs, features = "SM2_score11", group.by = "celltype", scale = TRUE) +coord_flip()

