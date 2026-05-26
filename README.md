# scRNA-Analysis-10X-Brain
This is an example of a preliminary analysis of a brain dataset from 10X genomics (https://www.10xgenomics.com/datasets/30-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-multiplexed-12-cm-os-3-1-standard-6-0-0)

## Quality Control
The first step to analyzing scRNA-seq data is to perform quality control. Here I look at the three metrics to determine which cells are low-quality and need to be removed from the dataset. 
- UMI Count - The total number of transcripts detected from a barcode
- Feature Count - The total nuber of features that a barcode has at least 1 transcript of
- Percent Mitochondrial Genes - Percentage of transcripts mapping to mitochondrial DNA

By making visualizations, I can check each QC metric and set thresholds to remove low quality cells. Here I will review each QC plot for the Brain_1 sample, and explain why we use them.

Log-log plot: This plot shows the UMI count for each barcode on the y-axis and orders all of the barcodes by their UMI count on the x-axis. In a good quality data set there will be a 'knee' separating the barcodes exposed to biological material and the empty ones.

<img width="400" height="400" alt="Screenshot 2026-05-26 at 10 48 47 AM" src="https://github.com/user-attachments/assets/6e3db67e-071c-4396-8413-9a91a0e64814" />

Histogram of UMI Count - This histogram visualizes the number of barcodes that have a certain number of total UMIs. Generally this will be a bimodal distribution, with a large group of empty barcodes below the UMI threshold and smaller group of good quality barcodes above the threhsold. We don't really see that here, but we can see that barcodes with UMI less than 1200 are clearly lower than most of the dataset.

<img width="400" height="400" alt="Screenshot 2026-05-26 at 10 55 00 AM" src="https://github.com/user-attachments/assets/f614f11a-edba-4dce-a143-fd6befdd1c1a" />

Histogram of Feature Count - This histogram typically has a very similar bimodal pattern to the total UMI count histogram described above. Here I set a threshold of 1000. There are just a couple hundred barcodes under this threshold, so I remove them to ensure we keep only good quality data.

<img width="400" height="400" alt="Screenshot 2026-05-26 at 11 00 27 AM" src="https://github.com/user-attachments/assets/eb14a583-c2d4-45ac-8e8c-114fdf5ad0f7" />

Histogram of Percent MT - In my experience, a typical dataset will have most of the data below 10-25% of mitochondrial transcripts. In order to not remove good cells, here I set the threshold at 25%.

<img width="400" height="400" alt="Screenshot 2026-05-26 at 11 03 58 AM" src="https://github.com/user-attachments/assets/147753f3-cf50-4315-8241-6e11372e57e7" />

Combination Scatterplot - This is typically the most helpful QC plot because I can visualize all three metrics onto one plot. We can clearly see the low-quality cells in the bottom left quadrant. 

<img width="400" height="400" alt="Screenshot 2026-05-25 at 6 48 50 PM" src="https://github.com/user-attachments/assets/17857b8a-814a-4250-aa51-4e11373ecc40" />

After removing these cells, I save the list of seurat objects and a table of final cells per sample.

## Clustering
After merging the seurat objects together I normalize, scale, and run PCA and UMAP. To decide if integration is needed, I look at the sample dimplots and look for evidence of batch effect. Below are the dimplots of each sample. There is good coverage across all three samples, so I move on to clustering the data. If there was evidence of batch effect correction, I would run Harmony to correct. 

<img width="400" height="400" alt="Screenshot 2026-05-25 at 6 56 51 PM" src="https://github.com/user-attachments/assets/c65f7445-8741-409d-aee6-8aceb5f6abe5" />

Now that the samples are combined and we addressed batch effect corretion, we cluster the data. The algorithm will look at all of the data and generate groups in an unsupervised manor, but I can adjust the resolution to generate more or less clusters. Here I set a resolution of 0.3.

<img width="400" height="400" alt="Screenshot 2026-05-25 at 6 59 18 PM" src="https://github.com/user-attachments/assets/b4bad3a3-dfef-4034-9e96-61ca0ed4890e" />

After I cluster the data, I perform differential expression to look for differences between the clusters. This heatmap shows the top ten genes by the average log2 fold-change for each cluster. There is pretty good distinction for most of the clusters. 

<img width="500" height="700" alt="Screenshot 2026-05-26 at 11 17 07 AM" src="https://github.com/user-attachments/assets/64fa485d-c480-40dc-8bb6-469939072756" />

In order to move forward to more downstream analyses, it is best to label the clusters. With the cluster QC plots and DE heatmap, the next step is to create marker gene plots. All three of these tools will be used to label the clusters.
