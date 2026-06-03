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
# Read in seurat object
ser <- readRDS('data/2_ser.RDS')

# Read in gene sets and convert to mouse
conv <- readRDS('/Users/michaelpatatanian/Desktop/Bioinformatics/gene_conv.RDS')
sets <- gmtPathways('data/h.all.v2026.1.Hs.symbols.gmt')
for(set_name in names(sets)){
  genes <- sets[[set_name]]
  mgi_genes <- conv$mgi[conv$hgnc %in% genes]
  mgi_genes <- mgi_genes[!is.na(mgi_genes) & mgi_genes != '']
  sets[[set_name]] <- mgi_genes
}
sets <- lapply(sets, unique)

for(clust in levels(ser$Cluster)){
  # Read in DE data
  marks <- read.csv(file.path(dir, paste0(clust, '_Brain_1.vs.Brain_2 DE.csv'))) %>%
    column_to_rownames('X')

  # Get ranks and run GSEA
  ranks <- -log(marks$p_val) * sign(marks$avg_log2FC)
  names(ranks) <- rownames(marks)
  ranks <- na.omit(ranks)
  gsea <- fgseaMultilevel(sets, ranks, eps=0)
  gsea <- gsea[order(-log(gsea$pval) * sign(gsea$NES), decreasing=T),]

  # Convert results to table and save
  max_le <- max(sapply(gsea$leadingEdge, length))
  mat <- matrix('', nrow=nrow(gsea), ncol=max_le + 6)
  rownames(mat) <- gsea$pathway
  colnames(mat) <- c(
    colnames(gsea)[-c(1,8)],
    paste0('Leading Edge ', 1:max_le))

  subset <- as.matrix(gsea[,2:7])
  mat[,1:6] <- subset
  for(r in 1:nrow(gsea)){
    le <- gsea$leadingEdge[r][[1]]
    le <- c(le, rep('', max_le - length(le)))
    mat[r, -c(1:6)] <- le
  }

  write.table(mat, file.path(dir, paste0(clust, '_Brain_1.vs.Brain_2 GSEA.csv')), 
    sep = ',', col.names=NA)
  
  # Make enrichment plots for top and bottom 5 gene sets
  pos <- head(gsea$pathway, 5)
  neg <- tail(gsea$pathway, 5)

  gsea_ylims <- c(min(gsea$ES)-0.1, max(gsea$ES))
  if(gsea_ylims[1] > 0) gsea_ylims[1] <- -gsea_ylims[2]
  if(gsea_ylims[2] < 0) gsea_ylims[2] <- -gsea_ylims[1]

  plots <- list()
  for(path in c(pos, neg)){
    # Grab gsea stats by gene
    stats <- rank(-ranks)
    ord <- order(stats)
    stats_adj <- ranks[ord]
    stats_adj <- sign(stats_adj) * (abs(stats_adj))
    stats_adj <- stats_adj/max(abs(stats_adj))
    pathway <- which(names(stats_adj) %in% sets[[path]])
    gsea_res <- calcGseaStat(stats_adj, selectedStats=pathway,
      returnAllExtremes=TRUE)
    
    # Convert to data frame x, y coordinates for plotting
    n <- length(stats_adj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(gsea_res$bottoms, gsea_res$tops))
    df <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

    # Generate plot
    p <- ggplot(df, aes(x=x, y=y)) + 
      geom_line(color='dark blue') + 
      geom_hline(yintercept=max(df$y), color='red', linetype='dashed') + 
      geom_hline(yintercept=min(df$y), color='red', linetype='dashed') + 
      geom_hline(yintercept=0, color='black') + 
      theme_bw() +
      geom_segment(mapping=aes(x=x, y=gsea_ylims[1], xend=x, yend=gsea_ylims[1]+0.05), size=0.2) +
      theme(panel.grid=element_blank()) +
      labs(x='Gene Rank', y='Enrichment Score', title=path)
  
    p <- p + ylim(gsea_ylims[1], gsea_ylims[2])

    # Add random permutations
    for(i in 1:100){
      rnk <- ranks
      names(rnk) <- sample(names(rnk))
      stats <- rank(-rnk)
      ord <- order(stats)
      stats_adj <- rnk[ord]
      stats_adj <- sign(stats_adj) * (abs(stats_adj))
      stats_adj <- stats_adj/max(abs(stats_adj))
      pathway <- which(names(stats_adj) %in% sets[[path]])
      gsea_res <- calcGseaStat(stats_adj, selectedStats=pathway,
        returnAllExtremes=TRUE)
      
      # Convert to data frame x, y coordinates for plotting
      n <- length(stats_adj)
      xs <- as.vector(rbind(pathway - 1, pathway))
      ys <- as.vector(rbind(gsea_res$bottoms, gsea_res$tops))
      tmp <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
      p$layers <- c(geom_line(data=tmp, color='light grey', size=.1), p$layers)
    }

    plots[[length(plots) + 1]] <- p
  }

  # Print plots
  pdf(file.path(dir, paste0(clust, '_Brain_1.vs.Brain_2_GSEA_enrichment.pdf')),
    height = 5, width = 5)
  for(plot in plots) print(plot)
  dev.off()
}