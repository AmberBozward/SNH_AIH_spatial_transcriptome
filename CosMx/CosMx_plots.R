# ==============================================
# Spatial Transcriptomics Figure Generation Script
# ==============================================
# Description: Generates UMAPs, barplots, boxplots, heatmaps, dotplots,
#              and CellChat analyses from a Seurat/CosMx object
# ==============================================

# ------------------------
# 0. Load required packages
# ------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(CellChat)
library(patchwork)
library(ggrepel)

# ------------------------
# 1. Load Seurat/CosMx object
# ------------------------
seurat_obj <- readRDS("data/seurat_obj_processed.rds")

# Create output folder if it doesn't exist
if(!dir.exists("figures")) dir.create("figures")

# ------------------------
# 2. UMAP Plots
# ------------------------
# UMAP by cell type
p1 <- DimPlot(seurat_obj, group.by = "final_cell_types") +
  ggtitle("UMAP by Cell Type")
ggsave("figures/UMAP_cell_types.png", p1, width=7, height=6)

# UMAP by sample
p2 <- DimPlot(seurat_obj, group.by = "sample") +
  ggtitle("UMAP by Sample")
ggsave("figures/UMAP_by_sample.png", p2, width=7, height=6)

# ------------------------
# 3. Barplots: Cell type proportions per sample
# ------------------------
cell_counts <- seurat_obj@meta.data %>%
  group_by(sample, final_cell_types) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count))

p_bar <- ggplot(cell_counts, aes(x=sample, y=prop, fill=final_cell_types)) +
  geom_bar(stat="identity") +
  theme_classic() +
  labs(y="Proportion", x="Sample", fill="Cell Type") +
  ggtitle("Cell Type Proportions per Sample")
ggsave("figures/barplot_cell_types.png", p_bar, width=8, height=6)

# ------------------------
# 4. Boxplots: Gene expression across cell types
# ------------------------
genes_to_plot <- c("CD63", "TIGIT", "IL17A") # modify as needed

for (gene in genes_to_plot) {
  p_box <- VlnPlot(seurat_obj, features = gene, group.by = "final_cell_types", pt.size=0.1) +
    ggtitle(paste0("Expression of ", gene))
  ggsave(paste0("figures/boxplot_", gene, ".png"), p_box, width=7, height=6)
}

# ------------------------
# 5. Heatmaps
# ------------------------
# Example: top 20 variable genes
top_genes <- head(VariableFeatures(seurat_obj), 20)
expr_matrix <- AverageExpression(seurat_obj, features = top_genes, group.by="final_cell_types")$RNA

pheatmap(expr_matrix, 
         filename = "figures/heatmap_top20.png", 
         cluster_rows=TRUE, cluster_cols=TRUE, 
         scale="row", fontsize_row=8)

# ------------------------
# 6. Dotplots
# ------------------------
p_dot <- DotPlot(seurat_obj, features = genes_to_plot, group.by = "final_cell_types") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle("Dotplot of Selected Genes")
ggsave("figures/dotplot_selected_genes.png", p_dot, width=8, height=6)

# ------------------------
# 7. CellChat Analysis
# ------------------------
# Only run if CellChat object has been set up previously
cellchat <- createCellChat(object = seurat_obj, group.by = "final_cell_types")
cellchat <- addMeta(cellchat, meta = seurat_obj@meta.data)
cellchat <- setIdent(cellchat, ident.use = "final_cell_types")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualize communication network (circle plot)
pdf("figures/cellchat_circle.pdf", width=8, height=8)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(seurat_obj$final_cell_types)),
                 weight.scale = T, label.edge= F, title.name = "CellChat Network")
dev.off()

# ==============================================
# End of Script
# ==============================================