Breast Cancer Single Cell RNAseq Analysis 

Abstract
Analyzing single-cell RNA-seq data from primary breast cancer tissue reveals the diversity of tumor and stromal cell populations. We focus on comparing the heterogeneity of Luminal A, Luminal B, HER2-positive, and Triple-negative breast cancer subtypes using this approach. Our study includes 14 patients and a total of 44,024 cells. We employ machine learning to cluster cells and identify robust marker genes distinguishing tumor from stromal cells. This process is repeated for each subtype, comparing highly expressed genes and calculating subtype fractions. Interestingly, we observed that heterogeneity related to cancer type is pronounced in tumor cells but not in stromal cell types.

Dataset
Our datasets originate from the Pan-cancer Blueprint datasets (http://blueprint.lambrechtslab.org). Among these, BC14 patient is classified as Luminal A, BC13 as Luminal B, BC01, BC06, BC08, and BC09 as HER2, while the remaining 8 patients are diagnosed with TNBC.

Methods
Single cell RNA sequence data preprocessing
We analyzed single-cell RNA sequence data from 44,024 cells isolated from primary breast cancer tissue, profiling a total of 33,694 genes using scanpy 1.6.1, a scalable toolkit for single-cell gene expression analysis. Quality control procedures were employed to filter out lowly expressed genes and low-quality cells. Specifically, genes detected in fewer than 3 cells were removed, and cells expressing fewer than 200 genes were filtered out. Additionally, cells with excessive expression of mitochondrial genes or total counts were excluded. Remarkably, our datasets exhibited high quality, with no cells filtered out, resulting in the retention of 26,040 genes.

Determination of major cell types 
We first normalized and log-transformed the data. Principal component analysis (PCA) was then conducted, followed by visualization using UMAP. Differential gene expression analysis was performed on clusters generated at various resolutions and using different numbers of neighbors and principal components (PCs) via t-tests and the Leiden algorithm. The top 5 highly expressed genes were identified for each cluster, aiding in the classification of cells into distinct cell types based on marker gene expression.

Subtyping for each cell type and phenotypic heterogeneity
Subtyping was carried out for each cell type using the same procedures as the previous step. Subgroups were delineated within each cell type based on marker gene expression, and the distribution of subgroups across different cancer types was analyzed. Heterogeneity within subgroups was then compared across the four breast cancer types.

Gene expression difference comparisons 
To identify differential gene expression across the four subtypes of breast cancer, highly expressed genes within each subgroup were identified based on previous clustering within each cell type. These genes were then selected and visualized in a heatmap to compare expression patterns among the four breast cancer types.

Results 
Please find results and figures from the results folder. 



