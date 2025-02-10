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

