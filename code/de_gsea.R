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

EnrichmentPlot <- function(stats, pathway, gsea_param, ylims=NULL, title=''){
  prep_df <- function(stats, pathway, gsea_param){
    # Get GSEA stats by gene
    rnk <- rank(-stats)
    ord <- order(rnk)
    stats_adj <- stats[ord]
    stats_adj <- sign(stats_adj) * (abs(stats_adj)^gsea_param)
    stats_adj <- stats_adj/max(abs(stats_adj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(stats_adj)))))
    pathway <- sort(pathway)
    gsea_res <- calcGseaStat(stats_adj, selectedStats=pathway,
      returnAllExtremes=TRUE)

    # Convert to data frame x, y coordinates for plotting
    bottoms <- gsea_res$bottoms
    tops <- gsea_res$tops
    n <- length(stats_adj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  }
  
  # Generate plot with real data
  df <- prep_df(stats, pathway, gsea_param)
  p <- ggplot(df, aes(x=x, y=y)) + 
    geom_line(color='dark blue') + 
    geom_hline(yintercept=max(df$y), colour='red', linetype='dashed') + 
    geom_hline(yintercept=min(df$y), colour='red', linetype='dashed') + 
    geom_hline(yintercept=0, colour="black") + 
    theme_bw() +
    geom_segment(mapping=aes(x=x, y=ylims[1], xend=x, yend=ylims[1]+0.05), size=0.2) +
    theme(panel.grid=element_blank()) +
    labs(x='Gene Rank', y='Enrichment Score', title=title)
  
  if(!is.null(ylims)) p <- p + ylim(ylims[1], ylims[2])


  # Add lines for random permutations
  for(i in 1:100){
    names(stats) <- sample(names(stats))
    df <- prep_df(stats, pathway, gsea_param)
    p$layers <- c(geom_line(data=df, color='light grey', size=.1), p$layers)
  }
  p
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