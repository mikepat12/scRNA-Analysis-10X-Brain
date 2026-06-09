# Run monocle3 container - podman run --rm -it monocle3:v1
# Add libraries
library(tidyverse)
library(Seurat)
library(monocle3)
if(!dir.exists('5 Pseudotime')) dir.create('5 Pseudotime')
dir <- '5 Pseudotime'

################################################################################
# Prep data
################################################################################
# Subset ser 
ser <- readRDS('data/2_ser.RDS')
ser <- ser[, !ser$Cluster %in% c('Pericyte', 'Granule', 'Macro', 'Vasc.Endo')]

# Create cds object
anno <- data.frame('gene_short_name' = rownames(ser))
rownames(anno) <- anno$gene_short_name
cds <- new_cell_data_set(ser[['RNA']]@counts, cell_metadata = ser[[]],
  gene_metadata = anno)

# Add umap
umap_coords <- Embeddings(ser, reduction = 'umap')
reducedDims(cds)$UMAP <- umap_coords

################################################################################
# Pseudotime
################################################################################
# Ensure factor and correct cell order
cluster <- factor(colData(cds)$Cluster)
names(cluster) <- colnames(cds)

# Single partition
partitions <- factor(rep(1, ncol(cds)))
names(partitions) <- colnames(cds)

# Build full cluster structure expected by monocle3
cds@clusters$UMAP <- list(clusters = cluster, partitions = partitions,
  resolution = NA, k = NA, method = "ser", num_clusters = length(levels(cluster))
)

# Learn trajectory 
cds <- learn_graph(cds, use_partition = TRUE)
saveRDS(cds, 'data/cds.RDS')
save_monocle_objects(cds, 'data/cds_monocle')

pdf(file.path(dir, '5 Trajectory.pdf'))
plot_cells(cds, color_cells_by = 'Cluster', label_groups_by_cluster=FALSE,
  label_leaves=FALSE, label_branch_points=FALSE)
dev.off()

# Order cells (Set neural progenitor cells as root)
root_cells <- colnames(cds)[colData(cds)$Cluster == 'NPC']
cds <- order_cells(cds, root_cells = root_cells)

# Plot pseudotime
pdf(file.path(dir, '5 Pseudotime.pdf'))
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups=FALSE, label_leaves=FALSE,
  label_branch_points=FALSE, graph_label_size=1.5)
dev.off()