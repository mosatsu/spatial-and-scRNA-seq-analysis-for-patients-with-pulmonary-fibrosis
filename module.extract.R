## Identify and extract the top genes with high S/N ratios that exhibit high expression levels within the SM2 group ##
modules = GetModules(seurat_obj)
library(WGCNA)
library(hdWGCNA)

# Extract all genes within the module
write.table(modules, file="module_gene.txt", sep="\t", quote=F, row.names=T, col.names=T)

# extract ...module_SM7
module_SM2 = modules[modules$module %in% "SM2",]
module_SM1 = modules[modules$module %in% "SM1",]
module_SM3 = modules[modules$module %in% "SM3",]
module_SM4 = modules[modules$module %in% "SM4",]
#....

# Extract the top 10 hub genes for each module
hub = GetHubGenes(seurat_obj, n_hubs=10)
write.table(hub, file="module_hubgene.txt", sep="\t", quote=F, row.names=T, col.names=T)

# Visualize the top 20 hub genes with a bubble plot
# Modify the annotation names to match the classification names
head(colnames(seurat_obj@meta.data))
p = DotPlot(seurat_obj, features=hub$gene_name, group.by="spot_category", dot.min=0.1) +
  coord_flip() + RotatedAxis() + scale_color_gradient2(high="red", mid="grey95", low="blue")
ggsave(p, file="hub_dotplot.png", dpi=300, width=6, height=15, device="png", limitsize=F, bg="white")

# Spatial visualization of genes
dir.create("features")

# Extract only the SM2 module
module_SM2 = modules[modules$module %in% "SM2",]
features = as.vector(module_SM2$gene_name)

for(i in c(1:nrow(module_SM2))){
  filename = paste0("./features/", features[i], ".png")
  p = SpatialFeaturePlot(seurat_obj, features=features[i])
  ggsave(p, file=filename, dpi=300, width=6, height=4, device="png", limitsize=F)
}



