# Add libraries
library(tidyverse)
library(Seurat)
library(Matrix)
library(scales)

################################################################################
# Helpers
################################################################################
read_data <- function(samp){
  library(Seurat)
  c <- file.path('/Users/michaelpatatanian/Desktop/cmc_test/', 
    samp, 'sample_feature_bc_matrix') %>%
    Seurat::Read10X()
  gx <- c[['Gene Expression']]
  cmo <- c[['Multiplexing Capture']]
  colnames(gx) <- paste0(colnames(gx), '_', samp)
  colnames(cmo) <- paste0(colnames(cmo), '_', samp)
  gx <- CreateSeuratObject(gx)
  gx[['CMO']] <- CreateAssayObject(cmo)
  return(gx)
}

make_breaks <- function(max){
  breaks <- seq(0, 2 * max, length.out=150)
  c(breaks, Inf)
}

################################################################################
# Quality Control
################################################################################
# Sample metadata
samples <- c('Brain_1', 'Brain_2', 'Brain_3', 'Brain_4')
metadata <- data.frame(Samples = samples, Tissue = 'Brain') %>%
  column_to_rownames('Samples')

# For each sample, read in data and create quality control plots
if(!dir.exists('1 Quality Control')) dir.create('1 Quality Control')
dir <- '1 Quality Control'
sers <- list()
samp_tab <- c()
for(samp in samples){
  # Read in data
  ser <- read_data(samp)

  # Calculate total UMI, features, and % mt genes
  ser <- PercentageFeatureSet(ser, '^mt-', col.name='Percent.MT')
  df <- data.frame('umi'=ser$nCount_RNA, 'nfeat'=ser$nFeature_RNA, 
    'pct.mt'=ser$Percent.MT, rank=rank(-ser$nCount_RNA))
  
  # Plot histograms of qc metrics
  plots <- list(
  # Log-log plot of UMI count
  ggplot(df, aes(x=rank, y=umi)) + 
    geom_line() +
    scale_x_continuous(trans='log10', 
      labels=label_number(suffix='k', scale=1e-3)) +
    scale_y_continuous(trans='log10', 
      labels=label_number(suffix='k', scale=1e-3)) +
    theme_bw() +
    xlab('Barcode Rank') +
    ylab('UMI Count') +
    geom_hline(yintercept=1200, color='red'),

  # Histogram of UMI count
  ggplot(df, aes(x=umi)) +
    geom_histogram(breaks=make_breaks(10000)) +
    theme_bw() +
    xlab('UMI Count') +
    ylab('# of Barcodes') +
    xlim(0, 10000) +
    labs(title='Histogram of UMI Count', subtitle=samp) +
    geom_vline(xintercept=1200, color='red'),

  # Histogram of Feature count
  ggplot(df, aes(x=nfeat)) +
    geom_histogram(breaks=make_breaks(4000)) +
    theme_bw() +
    xlab('Feature Count') +
    ylab('# of Barcodes') +
    xlim(0, 4000) +
    labs(title='Histogram of Feature Count', subtitle=samp) +
    geom_vline(xintercept=1000, color='red'),

  # Histogram of MT%
  ggplot(df, aes(x=pct.mt)) +
    geom_histogram(breaks=make_breaks(100)) +
    theme_bw() +
    xlab('% Mitochondrial Genes') +
    ylab('# of Barcodes') +
    xlim(0, 100) +
    labs(title='Histogram of % Mitochondrial Genes', subtitle=samp) +
    geom_vline(xintercept=25, color='red'),

  # QC scatter plot with all metrics
  ggplot(df, aes(x=umi, y=nfeat, color=pct.mt)) +
    geom_point(size=.2) +
    theme_classic() +
    xlim(0, 10000) +
    ylim(0, 4000) +
    geom_vline(xintercept=1200) +
    geom_hline(yintercept=1000) +
    xlab('UMI Count') +
    ylab('Feature Count') +
    labs(title='QC Scatter', subtitle=samp, color='% MT') +
    scale_color_gradient(low='grey', high='red') +
    theme(legend.position=c(0.2, 0.8))
  )

  # Print plots
  pdf(file.path(dir, paste0(samp, '-QC-plots.pdf')))
  for(plot in plots) print(plot)
  dev.off()

  # Filter counts based on thresholds
  filter <- df$umi > 1200 & 
    df$nfeat > 1000 & 
    df$pct.mt < 25
  
  ser <- ser[, filter]

  # Create metadata to add to ser
  cell_metadata <- metadata[rep(samp, ncol(ser)), , drop=F]
  rownames(cell_metadata) <- colnames(ser)
  ser@meta.data <- ser@meta.data %>%
    dplyr::select(-orig.ident) %>%
    mutate(Sample=samp) %>%
    cbind(cell_metadata)
  ser@meta.data <- dplyr::rename(ser@meta.data, UMI.Count=nCount_RNA,
    Feature.Count=nFeature_RNA)
  
  # Create table of sample cell numbers
  tmp_tab <- length(ser$Sample) %>%
    data.frame() %>%
    add_column(Sample = samp) %>%
    relocate(Sample) %>%
    dplyr::rename('Cell Numbers'='.')
  samp_tab <- rbind(samp_tab, tmp_tab)

  sers[[samp]] <- ser

}

# Save sample suerat objects and sample table
saveRDS(sers, 'data/1_sers.RDS')
write.table(samp_tab, file.path(dir, 'sample_cell_numbers.csv'), sep=',', row.names = F)