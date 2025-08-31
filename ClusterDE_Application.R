library(ClusterDE)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(Matrix)
library(hdf5r)
library(PairedData)
library(MASS)
library(scales)
library(readr)
library(magrittr)
library(knitr)
library(dplyr)

dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

#Random number generator, reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(2025)

## 1. data reading
h5_file <- "filtered_feature_bc_matrix_Bsub_minimal_media.h5"
message("Reading H5: ", h5_file)
mydata <- Read10X_h5(h5_file)
str(mydata)#genes by cells matrix
dim(mydata)#29765 genes & 2784 cells
dimnames(mydata)

# overview of data set
head(rownames(mydata), 20)#"aadK_1"  "aadK_2"  "aadK_3"  "aadK_4"
head(colnames(mydata), 20)#"AAACCTGAGCTTATCG-1" "AAACCTGAGGATGCGT-1"
as.matrix(mydata[1:8, 1:8])

# Create Seurat Object
obj <- CreateSeuratObject(counts = mydata)

## 2) Dataset overview (safe on sparse matrices) ----------------------------
count_mat <- GetAssayData(obj, slot = "counts")  # genes x cells (dgCMatrix)

# Basic dimensions
n_genes <- nrow(count_mat)
n_cells <- ncol(count_mat)

# Non-zero entries and sparsity
nnz <- length(count_mat@x)                         # number of nonzeros in sparse matrix
sparsity <- 1 - (nnz / (n_genes * n_cells))     # proportion of zeros: 98.61%

# Global UMI sum
total_umi <- sum(counts)

# Per-cell metrics (no densification)
umi_per_cell   <- Matrix::colSums(count_mat)
genes_per_cell <- Matrix::colSums(count_mat > 0)

# Per-gene metrics
umi_per_gene   <- Matrix::rowSums(counts)
cells_per_gene <- Matrix::rowSums(counts > 0)

# Light gene filter count (detected in ≥3 cells)
genes_ge3_cells <- sum(cells_per_gene >= 3)

# Helpful quantiles for narrative
cell_umi_q <- quantile(umi_per_cell,   probs = c(0.01, 0.25, 0.50, 0.75, 0.99), names = TRUE)
cell_feat_q<- quantile(genes_per_cell, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), names = TRUE)
gene_umi_q <- quantile(umi_per_gene,   probs = c(0.50, 0.75, 0.90, 0.95, 0.99), names = TRUE)
gene_cell_q<- quantile(cells_per_gene, probs = c(0.50, 0.75, 0.90, 0.95, 0.99), names = TRUE)

# Top 10 genes by total UMI
top10_idx <- order(umi_per_gene, decreasing = TRUE)[1:min(10, length(umi_per_gene))]
top10_tbl <- data.frame(
  Gene              = rownames(counts)[top10_idx],
  Total_UMI         = as.integer(umi_per_gene[top10_idx]),
  Detected_in_Cells = as.integer(cells_per_gene[top10_idx]),
  stringsAsFactors  = FALSE
)

## 3)  Table 1
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct <- function(x) sprintf("%.2f", 100 * x)

table1 <- data.frame(
  Metric = c(
    "Number of cells",
    "Number of genes (rows)",
    "Total UMI counts",
    "Non-zero entries",
    "Sparsity (%)",
    "Median UMI per cell",
    "Median detected genes per cell",
    "Genes detected in ≥3 cells"
  ),
  Value = c(
    fmt_int(n_cells),
    fmt_int(n_genes),
    fmt_int(total_umi),
    fmt_int(nnz),
    fmt_pct(sparsity),
    fmt_int(median(umi_per_cell)),
    fmt_int(median(genes_per_cell)),
    fmt_int(genes_ge3_cells)
  ),
  stringsAsFactors = FALSE
)

cell_summary <- data.frame(
  Quantile                  = names(cell_umi_q),
  UMI_per_cell              = as.integer(cell_umi_q),
  Detected_genes_per_cell   = as.integer(cell_feat_q),
  stringsAsFactors          = FALSE
)

gene_summary <- data.frame(
  Quantile               = names(gene_umi_q),
  Total_UMI_per_gene     = as.integer(gene_umi_q),
  Cells_detected_per_gene= as.integer(gene_cell_q),
  stringsAsFactors       = FALSE
)

# Save CSVs
write.csv(table1,       file = "results/table1_dataset_summary.csv", row.names = FALSE)
write.csv(cell_summary, file = "results/table1_cell_summary.csv",    row.names = FALSE)
write.csv(gene_summary, file = "results/table1_gene_summary.csv",    row.names = FALSE)
write.csv(top10_tbl,    file = "results/table1_top10_genes.csv",     row.names = FALSE)

# Markdown version 
md_lines <- c(
  "| Metric | Value |",
  "|:--|--:|",
  paste0("| ", table1$Metric, " | ", table1$Value, " |")
)
writeLines(md_lines, con = "results/table1_dataset_summary.md")


## Target Seurat pipeline + clustering + UMAP

RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj)
obj <- RunPCA(object = obj)
obj <- FindNeighbors(object = obj)
obj <- FindClusters(object = obj, resolution = 0.06)

# Export cluster sizes
cluster_sizes <- as.data.frame(table(Idents(obj)))
colnames(cluster_sizes) <- c("cluster", "n_cells")
write.csv(cluster_sizes, "results/target_cluster_sizes.csv", row.names = FALSE)

set.seed(2025)
obj <- RunUMAP(obj, dims = 1:10)

target_umap <- DimPlot(obj, reduction = "umap", label = TRUE) +
  NoLegend() + ggtitle("Target UMAP (resolution = 0.06)")
ggsave("results/umap_target_res006.png", target_umap, width = 6, height = 5, dpi = 300)

       
## Target DE analysis (0 vs 1) in Seurat
## - No thresholds: min.pct = 0, logfc.threshold = 0
## - Export: top list (Wilcoxon) + DE counts across 5 tests

#original_markers <- FindMarkers(obj,
#                                ident.1 = 0,
 #                               ident.2 = 1,
  #                              min.pct = 0,
   #                             logfc.threshold = 0) #
#print(paste("Original marker genes:", sum(original_markers$p_val_adj < 0.05)))

dir.create("results/target_DE", showWarnings = FALSE, recursive = TRUE)

## Baseline DE (Wilcoxon)
original_markers <- FindMarkers(
  object = obj,
  ident.1 = 0,
  ident.2 = 1,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0,
  
)

# Number of significant markers at BH 5%
n_sig_wilcox <- sum(original_markers$p_val_adj < 0.05, na.rm = TRUE)
message("Original marker genes (Wilcoxon, BH<0.05): ", n_sig_wilcox)

# Export a clean Top list (Top 20 by adjusted p-value)
top_wilcox <- original_markers
top_wilcox$Gene <- rownames(original_markers)
top_wilcox <- top_wilcox[order(top_wilcox$p_val_adj, decreasing = FALSE), 
                         c("Gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj")]
top_wilcox <- head(top_wilcox, 20)
write.csv(top_wilcox, "results/target_DE/table_top20_wilcox.csv", row.names = FALSE)

# Also export the full table for Wilcoxon
full_wilcox <- original_markers
full_wilcox$Gene <- rownames(original_markers)
full_wilcox <- full_wilcox[order(full_wilcox$p_val_adj, decreasing = FALSE), ]
write.csv(full_wilcox, "results/target_DE/DE_wilcox_full.csv", row.names = FALSE)

## 3.2 Compare 5 tests: wilcox / t / LR / poisson / negbinom
tests <- c("wilcox","t","LR","poisson","negbinom")
count_rows <- list()
all_results <- list()

for (te in tests) {
  mk <- FindMarkers(
    object = obj,
    ident.1 = 0,
    ident.2 = 1,
    test.use = te,
    min.pct = 0,
    logfc.threshold = 0,
    verbose = FALSE
  )
  mk$Gene <- rownames(mk)
  mk$test <- te
  mk <- mk[order(mk$p_val_adj, decreasing = FALSE), ]
  
  # Save per-test full table
  out_csv <- file.path("results/target_DE", paste0("DE_", te, "_full.csv"))
  write.csv(mk, out_csv, row.names = FALSE)
  
  # Count BH-significant genes
  n_sig <- sum(mk$p_val_adj < 0.05, na.rm = TRUE)
  count_rows[[length(count_rows) + 1]] <- data.frame(test = te, n_DE_BH_0_05 = n_sig)
  
  # Keep for optional overlap/plots
  all_results[[te]] <- mk
}

# Summary table: #DE by test
de_counts <- do.call(rbind, count_rows)
write.csv(de_counts, "results/target_DE/table_DE_counts_by_test.csv", row.names = FALSE)
print(de_counts)

# bar plot DE counts by test
de_counts$test <- factor(de_counts$test, levels = c("wilcox","t","poisson","negbinom","LR"))

p_bar <- ggplot(de_counts, aes(x = test, y = n_DE_BH_0_05, fill = test)) +
  geom_col(width = 0.7) +
  # centered labels on top of each bar
  geom_text(aes(label = comma(n_DE_BH_0_05)),
            vjust = -0.5, size = 4.2, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +         # harmonious, distinct colors
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +  # headroom for labels
  labs(title = "Number of DE genes across five tests (0 vs 1)",
       x = "Test", y = "# DE genes (BH < 0.05)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("results/target_DE/plot_DE_counts_by_test.png", p_bar,
       width = 7, height = 4.8, dpi = 300)
       
# table of "Top-20 DE genes (Wilcoxon)" 
in_csv  <- "results/target_DE/table_top20_wilcox.csv"
out_md  <- "results/target_DE/table_top20_wilcox.md"

tab <- read_csv(in_csv, show_col_types = FALSE)

stopifnot(all(c("Gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj") %in% names(tab)))

fmt_p <- function(x) {
  x <- ifelse(is.na(x), NA, x)
  # 处理极小p值
  out <- ifelse(x == 0, "<1e-300",
                ifelse(x < 1e-4, formatC(x, format = "e", digits = 2),
                       formatC(x, format = "f", digits = 4)))
  out
}

fmt_fc <- function(x) formatC(x, format = "f", digits = 2)
fmt_pct <- function(x) percent(x, accuracy = 0.1)  # 0.123 -> 12.3%

tab_fmt <- tab %>%
  mutate(
    `log2FC` = fmt_fc(avg_log2FC),
    `pct in 0` = fmt_pct(pct.1),
    `pct in 1` = fmt_pct(pct.2),
    `p-value`  = fmt_p(p_val),
    `q-value (BH)` = fmt_p(p_val_adj),
    Direction = ifelse(avg_log2FC > 0, "Up", "Down")
  ) %>%
  select(Gene, Direction, `log2FC`, `pct in 0`, `pct in 1`, `p-value`, `q-value (BH)`)

md_tbl <- kable(
  tab_fmt,
  format = "markdown",
  align = c("l","c","r","r","r","r","r"),
  caption = "Table 4.1. Top-20 DE genes (Wilcoxon) for cluster 0 vs cluster 1."
)
writeLines(md_tbl, out_md)

## Construct synthetic null data with ClusterDE

count_mat <- GetAssayData(object = obj, slot = "counts")
set.seed(2025)
system.time(
  synthetic_null <- constructNull(count_mat, 
                                  nCores = 1, 
                                  fastVersion = TRUE
                                  )
  )
# Save to file for reuse 
dir.create("results/null_data", showWarnings = FALSE, recursive = TRUE)
saveRDS(synthetic_null, "results/null_data/synthetic_null_matrix.rds")

# Quick sanity check
dim(synthetic_null)       # should match original matrix (genes × cells)
summary(as.vector(synthetic_null))[1:10]

## Synthetic null: Seurat pipeline + clustering + UMAP

#parallel the same pipeline 
RNGkind("L'Ecuyer-CMRG"); set.seed(2025)

#Create Seurat object from synthetic null matrix
obj_null <- CreateSeuratObject(counts = synthetic_null)

obj_null <- NormalizeData(object = obj_null)
obj_null <- FindVariableFeatures(object = obj_null)
obj_null <- ScaleData(object = obj_null)
obj_null <- RunPCA(object = obj_null)
obj_null <- FindNeighbors(object = obj_null)

#Clustering with resolution = 0.06
obj_null <- FindClusters(object = obj_null, resolution = 0.06)
obj_null <- RunUMAP(object = obj_null, dims = 1:10)
p1 <- DimPlot(object = obj_null, reduction = "umap", label = TRUE) + 
    ggtitle("Null UMAP (resolution = 0.06)") + NoLegend()

#Clustering with resolution = 0.075 (slightly higher granularity)
obj_null <- FindClusters(object = obj_null, resolution = 0.075)
obj_null <- RunUMAP(object = obj_null, dims = 1:10)
p2 <- DimPlot(object = obj_null, reduction = "umap", label = TRUE) + 
    ggtitle("Null UMAP (resolution = 0.075)") + NoLegend()

#Combine UMAP plots for comparison
p_combined <- p1+p2

ggsave(filename = "results/null_umap_comparison by resolution.png",
       plot = p_combined,
       width = 10, height = 5, dpi = 300)


#Null data DE analysis
null_markers <- FindMarkers(obj_null,
                            ident.1 = 0,
                            ident.2 = 1,
                            test.use = "wilcox",
                            min.pct = 0,        # same parameter
                            logfc.threshold = 0) # same

# Number of significant markers at BH 5%
n_sig_wilcox_null <- sum(null_markers$p_val_adj < 0.05, na.rm = TRUE)
message("null marker genes (Wilcoxon, BH<0.05): ", n_sig_wilcox_null)

#p_values and contrast scores
original_pval <- original_markers$p_val
names(original_pval) <- rownames(original_markers)

null_pval <- null_markers$p_val
names(null_pval) <- rownames(null_markers)
common <- intersect(names(original_pval), names(null_pval))

print(paste("Common genes for comparison:", length(common)))


res <- ClusterDE::callDE(original_pval[common], null_pval[common], 
                         nlogTrans = TRUE, FDR = 0.05)
cat(paste0("ClusterDE identified DE genes: ", length(res$DEgenes), "\n"))
cat(paste0("Original Seurat significant genes: ", sum(original_markers$p_val_adj < 0.05), "\n"))


#Correction by filter
original_markers_correct <- FindMarkers(
  object = obj,
  ident.1 = 0,
  ident.2 = 1,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  
)

# Number of significant markers at BH 5%
n_sig_wilcox_correct <- sum(original_markers_correct$p_val_adj < 0.05, na.rm = TRUE)
message("Original marker genes (Wilcoxon, BH<0.05, filtered): ", n_sig_wilcox_correct)

null_markers_correct <- FindMarkers(obj_null,
                            ident.1 = 0,
                            ident.2 = 1,
                            test.use = "wilcox",
                            min.pct = 0.1,        # same parameter
                            logfc.threshold = 0.25) # same

# Number of significant markers at BH 5%
n_sig_wilcox_null_correct <- sum(null_markers_correct$p_val_adj < 0.05, na.rm = TRUE)
message("null marker genes (Wilcoxon, BH<0.05,filtered): ", n_sig_wilcox_null_correct)

#p_values and contrast scores
original_pval_correct <- original_markers_correct$p_val
names(original_pval_correct) <- rownames(original_markers_correct)

null_pval_correct <- null_markers_correct$p_val
names(null_pval_correct) <- rownames(null_markers_correct)
common <- intersect(names(original_pval_correct), names(null_pval_correct))

print(paste("Common genes for comparison:", length(common)))


res_correct <- ClusterDE::callDE(original_pval_correct[common], null_pval_correct[common], 
                         nlogTrans = TRUE, FDR = 0.05)
cat(paste0("ClusterDE identified DE genes: ", length(res_correct$DEgenes), "\n"))
cat(paste0("Original Seurat significant genes: ", sum(original_markers_correct$p_val_adj < 0.05), "\n"))
res_correct$summaryTable

# Contrast score histogram (from res_correct)
p_cs <- ggplot(data = res_correct$summaryTable, aes(x = cs)) + geom_histogram(fill = "white", color = "black") + theme_bw() + ggtitle("Distribution of constrast scores")

# Save plot
dir.create("results/contrast_diag", showWarnings = FALSE, recursive = TRUE)
ggsave("results/contrast_diag/contrast_hist_filtered.png",
       plot = p_cs, width = 7, height = 5, dpi = 300)

## marker expression visualization
# Top 6 naive Seurat markers
p_naive <- FeaturePlot(
  obj[, obj$seurat_clusters %in% c(0, 1)], 
  features = rownames(original_markers_correct)[1:6], 
  ncol = 3
) + plot_annotation(title = "Naïve Seurat markers (Top 6)")

# Top 6 ClusterDE-calibrated markers
p_clusterde <- FeaturePlot(
  obj[, obj$seurat_clusters %in% c(0, 1)], 
  features = res_correct$summaryTable$Gene[1:6], 
  ncol = 3
) + plot_annotation(title = "ClusterDE-calibrated markers (Top 6)")

# Combine the two panels
p_combined <- p_naive / p_clusterde +
  plot_annotation(title = "Visualization of marker expression: Seurat vs ClusterDE")

#  Save to folder
dir.create("results/marker_viz", showWarnings = FALSE, recursive = TRUE)
ggsave("results/marker_viz/marker_expression_comparison.png",
       plot = p_combined, width = 10, height = 8, dpi = 300)


## Diagnostics: Symmetry of contrast scores

# Extract contrast scores from ClusterDE results
cs <- res_correct$summaryTable$cs
cs <- cs[is.finite(cs)]   # remove Inf/NA if any

# Split into positive and negative parts
pos <- cs[cs >= 0]
neg <- -cs[cs < 0]   # mirror negatives to positive scale

# If almost no negatives, this indicates asymmetry in itself
message("Number of positive cs: ", length(pos))
message("Number of negative cs: ", length(neg))

# Run Yuen's trimmed mean test if both sides have data
if (length(pos) > 5 & length(neg) > 5) {
  # Ensure equal length by random subsampling
  m <- min(length(pos), length(neg))
  set.seed(2025)
  pos_s <- sample(pos, m)
  neg_s <- sample(neg, m)
  
  # Apply Yuen’s t-test (20% trimming)
  yt <- yuen.t.test(pos_s, neg_s, trim = 0.2)
  print(yt)
} else {
  message("Not enough negative cs values for a balanced Yuen test.")
  yt <- NULL
}

sessionInfo()
