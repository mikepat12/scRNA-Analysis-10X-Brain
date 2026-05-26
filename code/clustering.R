# Add libraries
library(tidyverse)
library(Seurat)
library(harmony)
library(future)
set.seed(24)
options(future.globals.maxSize=Inf)

if(!dir.exists('2 Clustering')) dir.create('2 Clustering')
dir <- '2 Clustering'

# Read in seurat list and merge
sers <- readRDS('data/1_sers.RDS')
ser <- merge(sers[[1]], sers[-1])

# QC metrics to plot
qc_vars <- c('UMI.Count', 'Feature.Count', 'G2M.Score', 'S.Score', 
  'Percent.MT')

# Get genes for cell scoring
conv <- readRDS('/Users/michaelpatatanian/Desktop/Bioinformatics/gene_conv.RDS')
g2m <- conv$mgi[conv$hgnc %in% cc.genes$g2m.genes]
s <- conv$mgi[conv$hgnc %in% cc.genes$s.genes]

# Run processing without integration
plan(multicore, workers=3)
ser <- FindVariableFeatures(ser, verbose=T) %>%
  NormalizeData(verbose=T) %>%
  CellCycleScoring(s, g2m, verbose=T) %>%
  ScaleData(vars.to.regress=c('G2M.Score', 'S.Score', 'Percent.MT'), 
    verbose=T) %>%
  RunPCA(verbose=T) %>%
  RunUMAP(verbose=T, dims=1:30)
saveRDS(ser, 'data/2_pre_int_ser.RDS')

# Create sample dimplot
pdf(file.path(dir, 'sample_dimplots.pdf'))
DimPlot(ser, group.by = 'Sample', label = F, raster = F)
DimPlot(ser, group.by = 'Sample', split.by='Sample', ncol=2, raster=F) +
  theme(legend.position='none')
dev.off()

# QC plots without batch effect correction
pdf(file.path(dir, 'sample_qc_vars.pdf'))
for(qc in qc_vars){
  print(FeaturePlot(ser, qc, raster = F))
  print(VlnPlot(ser, qc, group.by = 'Sample', pt.size = 0) +
    theme(legend.position='none'))
}
dev.off()

# No batch effect correction needed

# Cluster data and plot clusters
ser <- FindNeighbors(ser, dims=1:30, verbose=F) %>%
  FindClusters(resolution=0.3, verbose=F)
ser$Cluster <- Idents(ser)

pdf(file.path(dir, 'clusters.pdf'))
DimPlot(ser, group.by = 'Cluster', label = T, raster = F)
dev.off()

# Violin plots of QC metrics by cluster
pdf(file.path(dir, 'cluster_qc_vars.pdf'), height = 4, width = 5)
for(qc in qc_vars){
  print(VlnPlot(ser, qc, group.by = 'Cluster', pt.size = 0) +
    theme(legend.position='none'))
}
dev.off()

# Run differential expression on broad clustering
marks <- FindAllMarkers(ser)
saveRDS(marks, 'data/all_cluster_marks.RDS')
write.table(marks, file.path(dir, 'All Cluster Markers.csv'), 
  sep=',', row.names=F, col.names=T)

# Make a heatmap of DE genes
top <- marks %>% group_by(cluster) %>% top_n(10, avg_log2FC)
width_factor <- case_when(length(unique(ser$Cluster)) < 3 ~ 4,
  length(unique(ser$Cluster)) < 3 ~ 4,
  length(unique(ser$Cluster)) < 5 ~ 2,
  TRUE ~ 1)
subser <- subset(ser, downsample=width_factor*200)

pdf(file.path(dir, 'de_heatmap.pdf'), width = ncol(subser) * width_factor / 300,
  height = nrow(top) / 10 + 1.5)
DoHeatmap(subser, top$gene) +
  theme(legend.position='none')
dev.off()

# Marker gene plots
marker_genes <- list(
  'Immune'=c('Ptprc'),
  'T Cells'=c('Cd3e', 'Cd4','Cd8a'),
  'NK'=c('Nkg7'),
  'B Cells'=c('Cd19', 'Ms4a1'),
  'Myeloid'=c('Cd14','Cd68','Itgam','Itgax','Csf3r'),
  'Endothelial'=c('Pecam1', 'Vwf', 'Lyve1'),
  'Epithelial'=c('Epcam','Cd24a'),
  'Pericyte'=c('Rgs5','Pdgfrb','Notch3'),
  'Fibroblasts'=c('Col1a1', 'Dcn'),
  'Neural'=c('Il17ra', 'Ikzf1', 'Calb1', 'Colq', 'Th', 'Avpr1a', 'Oprk1', 'Bmpr1b', 'Vcan', 'MrgA3', 'Asic3', 
    'Pvalb', 'Sst', 'Trpm8', 'Mrgprd', 'Mrgpra3', 'Ntrk1', 'Ntrk2', 'Ntrk3', 'Nefh', 
    'Calca', 'P2rx3', 'Cysltr2', 'Smr2')
)
if(!dir.exists('2 Clustering/marker_genes')) dir.create('2 Clustering/marker_genes')

for(i in 1:length(marker_genes)){
  # Grab name of list and get only genes in dataset
  name <- names(marker_genes)[i]
  gene_ids <- marker_genes[[i]] %in% rownames(ser)

  # For each list, create violin plots
  pdf(file.path(dir, paste0('marker_genes/', name, '_vln_markers.pdf')), height = 4, width = 5)
  for(mark in marker_genes[[i]][gene_ids]){
    print(VlnPlot(ser, mark, group.by = 'Cluster', pt.size = 0) +
      theme(legend.position = 'none'))
  }
  dev.off()
}

saveRDS(ser, 'data/2_ser.RDS')