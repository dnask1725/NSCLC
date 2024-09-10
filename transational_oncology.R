
library(biomaRt)
library(Seurat)
library(dplyr)
library(Matrix)

# Load your data
barcodes <- read.table("/GSM3635372_barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("/GSM3635372_genes.tsv", header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM("/GSM3635372_matrix.mtx")

# Set row names and column names
colnames(matrix) <- barcodes$V1
rownames(matrix) <- genes$V1

# Create Seurat object
sc_data <- CreateSeuratObject(counts = matrix, project = "NSCLC", min.cells = 3, min.features = 200)
print(sc_data)
# Map Ensembl IDs to gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_ensembl_ids <- rownames(sc_data)
gene_mapping_all <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                          filters = 'ensembl_gene_id', 
                          values = all_ensembl_ids, 
                          mart = ensembl)
View(gene_mapping_all)

# Create named vectors for mapping
ensembl_to_gene <- setNames(gene_mapping_all$hgnc_symbol, gene_mapping_all$ensembl_gene_id)
#gene_to_ensembl <- setNames(gene_mapping_all$ensembl_gene_id, gene_mapping_all$hgnc_symbol)

# Filter the Seurat object to keep only genes that have a mapping
mapped_genes <- rownames(sc_data) %in% names(ensembl_to_gene)
print(mapped_genes)
sc_data <- sc_data[mapped_genes, ]

View(sc_data)
library(openxlsx)
write.xlsx(gene_mapping_all, file = "/gene_mapping_all.xlsx", row.names = FALSE)
#__________________________________________________________________________________________________________________

# Create a named vector for mapping
ensembl_to_gene <- setNames(gene_mapping_all$hgnc_symbol, gene_mapping_all$ensembl_gene_id)
View(ensembl_to_gene)

# Replace Ensembl IDs with gene symbols in Seurat object
# Remove genes with no mapping

gene_names <- rownames(sc_data)
print(gene_names)
ensembl_ids <- names(ensembl_to_gene)
gene_symbols <- ensembl_to_gene
gene_mapping_df <- data.frame(Ensembl_ID = ensembl_ids, Gene_Symbol = gene_symbols, stringsAsFactors = FALSE)

View(gene_mapping_df)
# Save the dataframe to an Excel file
write.xlsx(gene_mapping_df, file = "/gene_mapping_processed.xlsx", row.names = FALSE)

#__________________________________________________
#FOR FILTERING ONLY ASSIGNED GENES

# Create a named vector for mapping

# Filter the Seurat object to keep only genes that have a mapping
mapped_genes <- rownames(sc_data) %in% names(ensembl_to_gene)
print(mapped_genes)
View(mapped_genes)
sc_data <- sc_data[mapped_genes, ]
gene_mapping_all <- gene_mapping_all[!is.na(gene_mapping_all$hgnc_symbol), ]
# Replace Ensembl IDs with gene symbols
gene_symbols <- ensembl_to_gene[rownames(sc_data)]
rownames(sc_data) <- gene_symbols

# Check for unassigned gene names (empty or NA)
unassigned_genes <- which(rownames(sc_data) == "" | is.na(rownames(sc_data)))
print(paste("Number of unassigned genes:", length(unassigned_genes)))

# Remove these genes from the Seurat object
if (length(unassigned_genes) > 0) {
  sc_data <- sc_data[-unassigned_genes, ]
}
print(paste("Dimensions after removing unassigned genes:", dim(sc_data)))

# Create a dataframe with Ensembl IDs and gene names of assigned genes
assigned_genes_df <- data.frame(
  Ensembl_ID = names(ensembl_to_gene)[match(rownames(sc_data), ensembl_to_gene)],
  Gene_Symbol = rownames(sc_data),
  stringsAsFactors = FALSE
)

# Save the dataframe to an Excel file
write.xlsx(assigned_genes_df, file = /assigned_genes.xlsx", row.names = FALSE)


#___________________________________________________

# 
# Quality Control and Filtering
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
print(sc_data)
VlnPlot(sc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
gene_names <- rownames(sc_data)
head(gene_names, 20)

# Filter cells based on QC metrics
# Adjust criteria based on the distribution of the metrics
sc_data <- CreateSeuratObject(counts = matrix, min.features = 200, min.cells = 3)

# Further subset cells based on the number of features
sc_data <- subset(sc_data, subset = nFeature_RNA > 200)
print(sc_data)
# Check dimensions to ensure cells are retained
print(dim(sc_data))

# Extract relevant QC data
qc_data <- data.frame(
  Cell_Barcode = colnames(sc_data),
  nFeature_RNA = sc_data$nFeature_RNA,
  nCount_RNA = sc_data$nCount_RNA,
  percent.mt = sc_data$percent.mt,
  stringsAsFactors = FALSE
)

# Save the QC data to an Excel file
write.xlsx(qc_data, file = "/qc_filtered_cells.xlsx", row.names = FALSE)
#__________________________________________________________


# Normalize the data
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
print(sc_data)
# Identify highly variable features
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
print(sc_data)
# Scale the data
all.genes <- rownames(sc_data)
sc_data <- ScaleData(sc_data, features = all.genes)

# Perform PCA
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))
print(sc_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc_data, dims = 1:2, reduction = "pca")
DimPlot(sc_data, reduction = "pca")

# Clustering the cells
sc_data <- FindNeighbors(sc_data, dims = 1:10)
sc_data <- FindClusters(sc_data, resolution = 0.5)
head(Idents(sc_data), 5)

# Run UMAP for visualization
sc_data <- RunUMAP(sc_data, dims = 1:10)
DimPlot(sc_data, reduction = "umap")

sc_data <- RunTSNE(sc_data, dims = 1:10)
DimPlot(sc_data, reduction = "tsne", label = TRUE)

# Find differentially expressed genes
# Find differentially expressed genes
sc_data.markers <- FindAllMarkers(sc_data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.01, test.use = "wilcox")
print(sc_data.markers)
write.xlsx(sc_data.markers, file = "/sc_markers.xlsx", row.names = FALSE)

markers <- nrow(significant_markers)
cat("Number of significant markers: ",markers, "\n")
significant_markers <- sc_data.markers[sc_data.markers$p_val_adj < 0.05, ]
print ()
# Get the top 10 genes for each cluster by average log2 fold change
top <- sc_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Print the top 10 genes for each cluster
print(top)
# Ensure the gene symbols are mapped to Ensembl IDs
# (This assumes you have already created the ensembl_to_gene named vector)
sc_data.markers$Ensembl_ID <- rownames(sc_data)[match(sc_data.markers$gene, rownames(sc_data))]
sc_data.markers$Gene_Symbol <- sc_data.markers$gene

# Select relevant columns (Ensembl IDs, gene symbols, and log2 fold change values)
marker_info <- sc_data.markers %>% select(Ensembl_ID, avg_log2FC)

# Save the dataframe to an Excel file
write.xlsx(marker_info, file = "/marker_info.xlsx", row.names = FALSE)


top <- sc_data.markers  %>% top_n(n = 10, wt = avg_log2FC)

min_pct <- 0.5  # Minimum 50% of cells in either cluster must express the gene
logfc_threshold <- 1  # Minimum log2 fold change of 1
# Function to find markers for a specific cluster
find_up_down_markers <- function(cluster_id, sc_data) {
  cluster_markers <- FindMarkers(sc_data, ident.1 = cluster_id, min.pct = min_pct, logfc.threshold = logfc_threshold, p_val_adj < 0.01)
  cluster_markers$gene <- rownames(cluster_markers)
  cluster_markers$cluster <- cluster_id
  return(cluster_markers)
}

all_clusters <- unique(Idents(sc_data))
all_markers <- lapply(all_clusters, find_up_down_markers, sc_data = sc_data)

# Combine the results into a single data frame
all_markers_df <- do.call(rbind, all_markers)

all_markers_df <- merge(all_markers_df, gene_mapping_all, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)


upregulated_genes <- sc_data.markers %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 2)

downregulated_genes <- sc_data.markers %>%
  filter(p_val_adj < 0.01 & avg_log2FC < -2)

# Print the number of up-regulated and down-regulated genes
cat("Number of up-regulated genes:", nrow(upregulated_genes), "\n")
cat("Number of down-regulated genes:", nrow(downregulated_genes), "\n")

library(biomaRt)

# Set up the Ensembl biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene annotation data, including gene type
annot <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
               mart = mart)

# View the annotation data
print(annot)
# Load the significant marker genes data from the saved Excel file
library(readxl)
file_path_upregu = "C:/Users/Saniya Kate/Desktop/PROJECTS/Seurat/upregulated.xlsx"
significant_markers_df = read_excel(file_path_upregu)

# Merge significant markers with annotation data
significant_markers_annot <- merge(significant_markers_df, annot, by = "ensembl_gene_id")
print(significant_markers_annot)
# Filter out pseudogenes
significant_markers_no_pseudogenes <- significant_markers_annot %>%
  filter(gene_biotype != "pseudogene")

# View the filtered significant markers
print(significant_markers_no_pseudogenes)

# Save the filtered significant markers to a new Excel file
write.xlsx(significant_markers_no_pseudogenes, file = "/upregulated_no_pseudogenes.xlsx", row.names = FALSE)

# Print the number of significant markers after removing pseudogenes
marker_genes <- nrow(significant_markers_no_pseudogenes)
cat("Number of significant markers after removing pseudogenes: ", marker_genes, "\n")


# Save the results
upregulated_genes <- merge(upregulated_genes, gene_mapping_all, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
View(upregulated_genes)
write.xlsx(upregulated_genes, file = "C:/Users/Saniya Kate/Desktop/PROJECTS/Seurat/upregulated.xlsx", row.names = FALSE)

downregulated_genes <- merge(downregulated_genes, gene_mapping_all, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
View(downregulated_genes)
write.xlsx(downregulated_genes, file = "/downregulated.xlsx", row.names = FALSE)



# Visualize top marker genes
top_genes <- bind_rows(upregulated_genes, downregulated_genes) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)
View(top_genes)
#________________________________________________________________


# Create a heatmap of top marker genes
DoHeatmap(sc_data, features = top$gene) + NoLegend()


#Violin
# Visualize top marker genes using a violin plot
top_genes <- unique(top$gene)

VlnPlot(sc_data, features = top_genes, ncol = 3)


# Save the top marker genes list
write.xlsx(top, file = "/top10_marker_genes.xlsx", row.names = FALSE)


# Save the Seurat object
saveRDS(sc_data, file = "sc_data.rds")


#________________________________________________________________

#GO and pathway analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("DOSE")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)


# Extract unique marker genes - so basically just top genes. 
marker_genes <- unique(top_genes$gene)
View(marker_genes)
# Map gene symbols to Entrez IDs
gene_list <- bitr(marker_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
View(gene_list)
#________________________________________________________________

# GO enrichment analysis
go_enrich <- enrichGO(gene = gene_list$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05)

# View results
head(go_enrich)
go_enrich_df <- as.data.frame(go_enrich)
go_enrich_df$geneID <- as.character(go_enrich_df$geneID)
go_enrich_df <- go_enrich_df %>%
  rowwise() %>%
  mutate(Gene_Symbols = paste(gene_list$SYMBOL[match(unlist(strsplit(geneID, "/")), gene_list$ENTREZID)], collapse = "/"))

# Save the merged results to an Excel file
write.xlsx(go_enrich_df, file = "/go_enrichment_results_with_gene_names.xlsx", row.names = FALSE)

#________________________________________________________________

# KEGG pathway analysis
kegg_enrich <- enrichKEGG(gene = gene_list$ENTREZID,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# View results
head(kegg_enrich)
kegg_enrich_df <- as.data.frame(kegg_enrich)
kegg_enrich_df$geneID <- as.character(kegg_enrich_df$geneID)
kegg_enrich_df<- kegg_enrich_df %>%
  rowwise() %>%
  mutate(Gene_Symbols = paste(gene_list$SYMBOL[match(unlist(strsplit(geneID, "/")), gene_list$ENTREZID)], collapse = "/"))

# Save the results to an Excel file
write.xlsx(kegg_enrich_df, file = "/kegg_enrichment_results.xlsx", row.names = FALSE)

#________________________________________________________________
#Vizualizing GO and KEGG:
# Load necessary libraries
library(clusterProfiler)
library(enrichplot)

# Dotplot for GO enrichment
dotplot(go_enrich, showCategory = 10)

# Barplot for KEGG enrichment
barplot(kegg_enrich, showCategory = 30)

# Visualize the top 10 GO terms using a bar plot
barplot(go_enrich, showCategory = 10, title = "Top 10 GO Enrichment Terms")

# Visualize the top 10 GO terms using a dot plot
dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Enrichment Terms")

#________________________________________________________________________________
#PPI NETWORK
library(STRINGdb)
library(dplyr)
library(openxlsx)

# Read upregulated and downregulated genes from files
upregulated_genes <- read.xlsx("/upregulated.xlsx")
View(upregulated_genes)
downregulated_genes <- read.xlsx("/downregulated.xlsx")

# Add status column
upregulated_genes$status <- "upregulated"
downregulated_genes$status <- "downregulated"

upregulated_genes <- upregulated_genes[, c("gene", "status", "ensembl_gene_id")]
View(upregulated_genes)
downregulated_genes <- downregulated_genes[, c("gene", "status", "ensembl_gene_id")]

degs <- rbind(upregulated_genes, downregulated_genes)

# Rename the column
names(degs)[names(degs) == "gene"] <- "GeneSymbol"

# Ensure the dataframe `degs` is correctly structured
str(degs)

ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)

# Convert to data frames if needed
ppi_network <- as.data.frame(ppi_network)
mapped_genes <- as.data.frame(mapped_genes)

# Merge PPI network with gene names for 'from' and 'to' columns
ppi_network <- merge(ppi_network, mapped_genes[, c("STRING_id", "GeneSymbol")], by.x = "from", by.y = "STRING_id", all.x = TRUE)
names(ppi_network)[names(ppi_network) == "GeneSymbol"] <- "from_gene"
ppi_network <- merge(ppi_network, mapped_genes[, c("STRING_id", "GeneSymbol")], by.x = "to", by.y = "STRING_id", all.x = TRUE)
names(ppi_network)[names(ppi_network) == "GeneSymbol"] <- "to_gene"

ppi_network <- ppi_network[, c("from_gene", "to_gene", "from", "to", "combined_score")]
head(ppi_network)

write.table(ppi_network, file = "ppi_network_with_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


#____________________________________________________________________________
# CYTOSCAPE
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RCy3")

# Load the RCy3 library
library(RCy3)

cytoscapePing()
ppi_network_file <- "/ppi_network_with_genes.txt"

# Check if the file exists
if (!file.exists(ppi_network_file)) {
  stop("File not found: ", ppi_network_file)
}

# Start Cytoscape and connect
cytoscapePing()

# Import the table into Cytoscape
tryCatch({
  res <- commandsPOST('table from file', list(file = ppi_network_file, firstRowAsColumnNames = TRUE))
  cat("Table imported successfully\n")
}, error = function(e) {
  cat("Error importing table: ", e$message, "\n")
})

# Create the network from the table by specifying source and target columns
tryCatch({
  res <- commandsPOST('network create from table', list(
    sourceColumn = "from_gene",
    targetColumn = "to_gene",
    interactionColumn = "combined_score",
    sourceTable = "node:default",
    targetTable = "edge:default"
  ))
  cat("Network created successfully\n")
}, error = function(e) {
  cat("Error creating network: ", e$message, "\n")
})
