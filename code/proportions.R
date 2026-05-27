# Add libraries
library(tidyverse)
library(Seurat)

if(!dir.exists('3 Proportions')) dir.create('3 Proportions')
dir <- '3 Proportions'

# Read in seurat object
ser <- readRDS('data/2_ser.RDS')

################################################################################
# Proportions
################################################################################
# Calculate clusters as proportion of each sample
pdf(file.path(dir, 'cluster_proportions.pdf'), height = 5, width = 5)
for(clust in Idents(ser$Cluster)){
  # Get proportions of cluster by sample
  ss <- subset(ser, idents=clust)
  tab <- table(ss$Sample)
  tab <- signif(100*tab/table(ser$Sample), 3)

  # Create colplot of samples
  df <- data.frame(tab)
  colnames(df) <- c('Sample', 'Proportion')
  p <- ggplot(df, aes(x=Sample, y=Proportion, fill=Sample)) +
    geom_col() +
    theme_bw() +
    xlab('') +
    ylab('Proportion of total (%)') +
    labs(title=paste('Cluster', clust)) +
    scale_fill_brewer(palette='Set1') +
    theme(legend.position='none',
      axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
  print(p)
}
dev.off()

# Save table of cluster by sample
table(ser$Sample, ser$Cluster) %>%
  write.table(file.path(dir, 'Cluster-Sample Cell Numbers.csv'), 
    sep=',', col.names=NA)