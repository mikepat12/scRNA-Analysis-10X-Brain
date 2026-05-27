# Add libraries
library(tidyverse)
library(Seurat)
library(EnhancedVolcano)
library(fgsea)
library(future)
options(future.globals.maxSize=Inf)
plan(multicore, workers = 2)

if(!dir.exists('4 DE GSEA')) dir.create('4 DE GSEA')
dir <- '4 DE GSEA'

################################################################################
# Helpers
################################################################################
# For volcano plot axes
get_min_max <- function(vals){
  vals <- sort(vals)
  bot_diffs <- sapply(2:11, function(x) vals[x] - vals[x-1])
  bot_mids <- which(c(T, bot_diffs > 10 * sd(vals)))
  bot_thresh <- vals[max(bot_mids)] - sd(vals)

  vals <- sort(vals, decreasing=T)
  top_diffs <- sapply(2:11, function(x) vals[x] - vals[x-1])
  top_mids <- which(c(T, -top_diffs > 10 * sd(vals)))
  top_thresh <- vals[max(top_mids)] + sd(vals)
  return(c(bot_thresh, top_thresh))
}

################################################################################
# Differential Expression and GSEA
################################################################################
# Read in seurat object
ser <- readRDS('data/2_ser.RDS')

# Compare Brain_1 and Brain_2 in each cluster
for(clust in levels(ser$Cluster)){
  # Grab cells to compare
  cells_1 <- colnames(ser)[ser$Cluster == clust & ser$Sample == 'Brain_1']
  cells_2 <- colnames(ser)[ser$Cluster == clust & ser$Sample == 'Brain_2']

  # Calculate differential expression
  marks <- FindMarkers(ser@assays[['RNA']], cells.1=cells_1, cells.2=cells_2, 
    logfc.threshold=0, features=NULL, verbose=F, test.use='wilcox')

  # Add pvalue correction and save results
  marks$p_val[marks$p_val == 0] <- 1e-310
  marks$p_val_adj <- p.adjust(marks$p_val, method='BH')
  write.table(marks, file.path(dir, paste0(clust, '_Brain_1.vs.Brain_2 DE.csv')), 
    sep=',', col.names=NA)

  # Create volcano plot
  pdf(file.path(dir, paste0(clust, '_Brain_1.vs.Brain_2_volcano.pdf')),
    height = 5, width = 5)
  res <- marks[!is.na(marks$p_val_adj) & !is.na(marks$avg_log2FC), ]
  res$p_val_adj[res$p_val_adj == 0] <- 1e-310
  x_lims <- get_min_max(res$avg_log2FC)
  y_lims <- get_min_max(-log10(res$p_val_adj))
  if(y_lims[1] < 0) y_lims[1] <- 0
  if(y_lims[2] > 315) y_lims[2] <- 315

  p <- EnhancedVolcano(res, rownames(res), 'avg_log2FC', 'p_val_adj', 
    FCcutoff=0, , xlab=bquote(~Log[2] ~ "average fold change"), 
      ylab=bquote(~-Log[10] ~ italic(p_val_adj)), title=paste('Cluster', clust),
      legendPosition='none') +
    annotate(geom='text', label='Brain_1', y=Inf, x=Inf, size=5,
      hjust=1, vjust=1) +
    annotate(geom='text', label='Brain_2', y=Inf, x=-Inf, size=5, 
      hjust=0, vjust=1) +
    xlim(x_lims[1], x_lims[2]) +
    ylim(y_lims[1], y_lims[2]) +
    theme_bw() +
    theme(legend.position='none', panel.grid=element_blank()) +
    geom_hline(yintercept=0)
  print(p)
  dev.off()
}

################################################################################
# Gene set enrichment analysis
################################################################################