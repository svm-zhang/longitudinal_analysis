library(Seurat)
library(data.table)
library(ggplot2)

plot_theme <- theme(
  axis.line = element_line(color = "black"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
)

# set up folder where plots are dumped.
plot_dir <- "./sc_plots"
if (! dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# function to identify DEGs between specified clusters.
find_degs_between_clusters <- function(obj, c1, c2) {
  dt <- as.data.table(
    FindMarkers(obj, ident.1 = c1, ident.2 = c2, only.pos=TRUE),
    keep.rownames=TRUE
  )
  names(dt)[1] <- "gene"
  # I would also look at fold change in addition to p value in real practice.
  dt[p_val_adj <= 0.05]
}

# Load the data
data_folder <- "./filtered_gene_bc_matrices/hg19/"
data <- Read10X(data.dir = data_folder)

pbmc <- CreateSeuratObject(counts = data, min.cells=3)

# get percent of mt genes per cell
pbmc[["PctMITO"]] <- PercentageFeatureSet(pbmc, pattern="^MT-")
# This metric tries to capture the complexity per cell.
pbmc[["GenePerUMI"]] <- log(pbmc$nFeature_RNA) / log(pbmc$nCount_RNA)
# pbmc@metedata
# Each row represents one cell, row name is barcode ID.
# columns: nCount_RNA, nFeature_RNA, pct_mt
metadata_dt <- pbmc@meta.data

# Visualize the distribution of UMI counts.
# VlnPlot() function from seurat can also serve the purpose but I like
# simple histogram.
# Majority cells in the data set have at least UMI count of 1000.
# The cell has the least UMI captureed is 546 in this dataset.
m <- ggplot(metadata_dt, aes(x=nCount_RNA)) +
  geom_histogram(binwidth = 100) +
  geom_vline(xintercept = 500, linetype = "dashed", color="red") +
  geom_vline(xintercept = 1000, linetype = "dashed", color="red") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "data.umi.hist.pdf"), m, height=10, width=14)

# Visualize the distribution of uniquely captured genes.
m <- ggplot(metadata_dt, aes(x=nFeature_RNA)) +
  geom_histogram(binwidth = 100) +
  geom_vline(xintercept = 200, linetype = "dashed", color="red") +
  geom_vline(xintercept = 2000, linetype = "dashed", color="red") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "data.gene.hist.pdf"), m, height=10, width=14)

# Visualize the distribution of percentage of captured mitochondria genes.
m <- ggplot(metadata_dt, aes(x=PctMITO)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 6, linetype = "dashed", color="red") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "data.pct_mito.hist.pdf"), m, height=10, width=14)

# If GenePerUMI is low, it means low complexity. In other words, we capture
# a low number of genes and simply sequenced these limited number of genes
# over and over again. Cell such as red blood cell or dying cells can have
# such profile. Usually I find 0.8 is a good threshold.
# And this specific dataset does not seem to have concerning complexity issue.
m <- ggplot(metadata_dt, aes(x=GenePerUMI)) +
  geom_histogram(binwidth = 0.001) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color="red") +
  xlim(0.5, 1) +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "data.gene_per_umi.hist.pdf"), m, height=10, width=14)

# Visualize relationship between UMI and gene counts colored by percentage of
# mitochondria genes. thresholds of 200-2000 for nFeature_RNA and 500 for
# nCount_RNA seem to be good ones.
# FeatureScatter() function from seurat serves similar purpose.
# Potiential problems to look for:
# bottom left quadrant: cells with low number of UMI and genes
# bottom right quadrant: low complexity cells
# slope: we expect a strong positive correlation
m <-ggplot(metadata_dt, aes(x=nCount_RNA, y=nFeature_RNA, color=PctMITO)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500, linetype = "dashed", color="red") +
  geom_hline(yintercept = 200, linetype = "dashed", color="red") +
  geom_hline(yintercept = 2000, linetype = "dashed", color="red") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "data.umi_vs_gene.pdf"), m, height=10, width=14)

# filter original dataset by criteria:
# 1. cells with low UMI
# 2. cells with low number of unique gene captured
# 3. cells with high percentage of mitochondria gene captured
# 4. cells with low complexity
# Sometimes, if we see abnormally large value of nCount_RNA, we also
# need to check possibility of doublets and/or multiplets.
pbmc <- subset(
  pbmc, subset = (nCount_RNA >= 500) &
    (nFeature_RNA >= 200) &
    (nFeature_RNA <= 2000) &
    (PctMITO < 6) &
    (GenePerUMI >= 0.8)
)
saveRDS(file = "q4.rds", pbmc)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# determine number of components
m <- ElbowPlot(pbmc)
ggsave(file.path(plot_dir, "pca.elbow.pdf"), m, height=12, width=16)
n_components <- 9

pbmc <- FindNeighbors(pbmc, dims = 1:n_components)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:n_components)
m <- DimPlot(pbmc, reduction = "umap")
ggsave(file.path(plot_dir, "cell.cluster.umap.pdf"), m, height=12, width=16)

# identify DEGs that differ between cluster 8 and 2. 
deg_res_clusters_dt <- find_degs_between_clusters(pbmc, 8, 2)

# if the purpose is to identify DEGs between one cell cluster versus all
# remaining cluster, use the following.
deg_one_vs_all_dt <- as.data.table(
  FindAllMarkers(pbmc, only.pos = TRUE), keep.rownames=TRUE
)
names(deg_one_vs_all_dt)[1] <- "gene"
# names(deg_one_vs_all_dt)
# p_val, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, cluster, gene
deg_one_vs_all_dt <- deg_one_vs_all_dt[p_val_adj <= 0.05]

