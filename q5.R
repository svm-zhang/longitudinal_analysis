library(Seurat)
library(SeuratData)
library(ggplot2)
library(data.table)

InstallData("ifnb")

plot_theme <- theme(
  axis.line = element_line(color = "black"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
)
plot_dir <- "./q5_plots"
if (! dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

ifnb <- LoadData("ifnb")
ifnb[["PctMITO"]] <- PercentageFeatureSet(ifnb, pattern="^MT-")
ifnb[["GenePerUMI"]] <- log(ifnb$nFeature_RNA) / log(ifnb$nCount_RNA)
metadata_dt <- ifnb@meta.data

# I am not including the plot code that guided me with the filter decision.
# Please refer to the q4.R.
# I do include the pdfs for these plots in the expected_output folder.
# Apply filters. Percent of mitochondria gene is irrelevant in this data.
ifnb <- subset(ifnb,, subset = (nCount_RNA >= 500) &
  (nFeature_RNA >= 200) &
  (nFeature_RNA <= 2000) &
  (GenePerUMI >= 0.8)
)
# split expression measurement by condition: control versus interferon stimulation
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst")
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
# review the elbow plot and determine number of components
m <- ElbowPlot(ifnb)
ggsave(file.path(plot_dir, "pca.elbow.pdf"), m, height=12, width=16)
n_components <- 15
# cluster cells regardless of condition, a.k.a no integration
ifnb <- FindNeighbors(ifnb, dims = 1:n_components)
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:n_components, reduction.name = "umap.unintegrated")

# integrate same cell types from the 2 conditions
print("Integrate cell types from control and interferon stimulated groups.")
ifnb <- IntegrateLayers(
  object = ifnb,
  method = CCAIntegration, orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

# cluster cells after alignment
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:n_components)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:n_components, reduction = "integrated.cca")

# prepare cluster data and embeddings for visualization
cluster_with_embed_dt <- as.data.table(
  ifnb@meta.data
)[, c("stim", "seurat_annotations", "seurat_clusters")]
cluster_with_embed_dt[, cell_id := rownames(ifnb@meta.data)]
setkey(cluster_with_embed_dt, cell_id)
# grab embedding before and after integration
embed_unintg <- Embeddings(ifnb, reduction = "umap.unintegrated")
embed_intg <- Embeddings(ifnb, reduction = "umap")
embed_dt <- as.data.table(cbind(embed_unintg, embed_intg))
embed_dt[, cell_id := rownames(embed_unintg)]
setkey(embed_dt, cell_id)
cluster_with_embed_dt <- embed_dt[cluster_with_embed_dt]

# visualize cell cluster assignment for different condition before and after
# integration
m <- ggplot(cluster_with_embed_dt, aes(x=umapunintegrated_1, y=umapunintegrated_2, color=stim)) +
  geom_point(size=0.5) +
  theme_bw() + plot_theme
ggsave(
  file.path(plot_dir, "cell.cluster.by_cond.prior_integration.pdf"),
  m,
  height=12,
  width=16
)
m <- ggplot(cluster_with_embed_dt, aes(x=umap_1, y=umap_2, color=stim)) +
  geom_point(size=0.5) +
  theme_bw() + plot_theme
ggsave(
  file.path(plot_dir, "cell.cluster.by_cond.post_integration.pdf"),
  m,
  height=12,
  width=16
)
m <- ggplot(cluster_with_embed_dt, aes(
  x=umapunintegrated_1, y=umapunintegrated_2, color=seurat_annotations
  )) +
  geom_point(size=0.5) +
  facet_wrap(~ stim, ncol=2) +
  theme_bw() + plot_theme
ggsave(
  file.path(plot_dir, "cell.cluster.by_type.prior_integration.pdf"),
  m,
  height=12,
  width=16
)
m <- ggplot(cluster_with_embed_dt, aes(x=umap_1, y=umap_2, color=seurat_annotations)) +
  geom_point(size=0.5) +
  facet_wrap(~ stim, ncol=2) +
  theme_bw() + plot_theme
ggsave(
  file.path(plot_dir, "cell.cluster.by_type.post_integration.pdf"),
  m,
  height=12,
  width=16
)

# now we identify DEGs between condition for each cell type
Idents(ifnb) <- paste(ifnb$seurat_annotations, ifnb$stim, sep="_")
cell_types <- unique(as.character(ifnb$seurat_annotations))
degs_by_cond_per_cell <- list()
# this takes time per cell type.
# break earlier if needed.
for (i in seq_len(length(cell_types))){
  cell_type <- cell_types[i]
  print(
    paste("Compare DEGs between condition for cell type ", cell_type, sep="")
  )
  deg_dt <- as.data.table(FindMarkers(
      ifnb,
      ident.1 = paste(cell_type, "STIM", sep="_"),
      ident.2 = paste(cell_type, "CTRL", sep="_"),
      verbose=FALSE
  ), keep.rownames = TRUE)
  names(deg_dt)[1] <- "gene"
  deg_dt <- deg_dt[p_val_adj <= 0.05]
  deg_dt[, cell_type := cell_type]
  print(head(deg_dt))
  degs_by_cond_per_cell[[cell_type]] <- deg_dt
}
saveRDS(file = "q5.deg.rds", degs_by_cond_per_cell)

