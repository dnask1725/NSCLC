# NSCLC

## Data Colection
The publicly available scRNA-seq data were downloaded from the GEO database. Data processing was performed using the Seurat package in R. For quality control check, I removed genes with a minimum number of features 200 having non zero counts and a minimum number of cells as 3. The filtered data was normalized using log-transformation and was used further.

## Differental Gene Expression
I identified DEGs using the Wilcoxon rank sum test,  based on the threshold adj.P.Val < 0.01 and logFC > 2 for up-regulated genes and adj.P.Val < 0.01 and logFC < -2 for down-regulated genes.

## GO-KEGG enrichment
The GO functional and KEGG pathway enrichment analyses associated with the DEGs was performed using analysis tools


## Results
