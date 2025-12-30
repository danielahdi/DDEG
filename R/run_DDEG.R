#' DDEG: Differential Data Expression Generator
#' @importFrom stats as.formula reorder
#' @importFrom utils head read.csv read.table write.csv
#' @param count_file Path to the count matrix (CSV or TXT)
#' @param meta_file (Optional) Path to the metadata file
#' @param group_column (Optional) The column name in metadata
#' @param p_val_cutoff Cutoff for adjusted p-value (default 0.05)
#' @param logfc_cutoff Cutoff for absolute log2 fold change (default 1)
#' @export
run_DDEG <- function(count_file, meta_file = NULL, group_column = NULL, p_val_cutoff = 0.05, logfc_cutoff = 1) {
  message(">>> Loading libraries and data...")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message(">>> Installing gene annotation database...")
    BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  }
  # 1. Read count data
  read_data <- function(path) {
    if (grepl("\\.csv", path)) return(read.csv(path, row.names = 1, check.names = FALSE))
    return(read.table(path, header = TRUE, row.names = 1, check.names = FALSE))
  }
  counts <- read_data(count_file)
  # 2. Pre-filtering: Keep genes with at least 10 counts in total
  message(">>> Pre-filtering: Removing genes with total counts < 10...")
  counts <- counts[rowSums(counts) >= 10, ]
  # 3. Clean Gene IDs and round decimals
  message(">>> Cleaning Gene IDs and preparing counts...")
  rownames(counts) <- gsub("\\..*", "", rownames(counts))
  counts <- round(counts)
  # 4. Handle Metadata
  if (is.null(meta_file)) {
    col_names <- colnames(counts)
    group_labels <- rep("Treatment", length(col_names))
    control_idx <- grep("control|ctrl|normal|healthy|wild|wt", col_names, ignore.case = TRUE)
    if (length(control_idx) == 0) {
      half <- ceiling(length(col_names) / 2)
      group_labels[1:half] <- "Group1"; group_labels[(half+1):length(col_names)] <- "Group2"
    } else { group_labels[control_idx] <- "Control" }
    metadata <- data.frame(condition = factor(group_labels), row.names = col_names)
    target_col <- "condition"
  } else {
    metadata <- read_data(meta_file)
    target_col <- if (is.null(group_column)) colnames(metadata)[1] else group_column
    common_samples <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, common_samples]
    metadata <- metadata[common_samples, ]
  }
  # 5. Run DESeq2 Analysis
  message(paste(">>> Running DESeq2 analysis grouped by:", target_col))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                        design = as.formula(paste0("~", target_col)))
  dds <- DESeq2::DESeq(dds)
  # تغییر: خاموش کردن فیلترینگ مستقل و تشخیص outliers برای کاهش NA
  res <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
  deg_results <- as.data.frame(res)
  # 6. Adding Gene Symbols
  message(">>> Mapping Ensembl IDs to Gene Symbols...")
  deg_results$gene_symbol <- tryCatch({
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                          keys = rownames(deg_results),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  }, error = function(e) {
    rep(NA, nrow(deg_results))  # در صورت خطا، همه NA
  })
  # جدید: جایگزینی NA در gene_symbol با Ensembl ID
  deg_results$gene_symbol[is.na(deg_results$gene_symbol)] <- rownames(deg_results)[is.na(deg_results$gene_symbol)]
  # 7. NA Management and Status Setting
  # جایگزینی NAها با مقادیر مناسب برای تمیز کردن خروجی
  deg_results$log2FoldChange[is.na(deg_results$log2FoldChange)] <- 0
  deg_results$pvalue[is.na(deg_results$pvalue)] <- 1
  deg_results$padj[is.na(deg_results$padj)] <- 1
  deg_results$diff <- "Stable"
  deg_results$diff[deg_results$log2FoldChange > logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Up"
  deg_results$diff[deg_results$log2FoldChange < -logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Down"
  # 8. Output Preparation & Saving
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- paste0("DDEG_Results_", timestamp)
  if (!dir.exists(output_dir)) dir.create(output_dir)
  deg_results$ensembl_id <- rownames(deg_results)
  final_cols <- c("ensembl_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "diff")
  deg_results_output <- deg_results[, final_cols]
  write.csv(deg_results_output, file.path(output_dir, "DEG_Full_Table.csv"), row.names = FALSE)
  # 9. Visualization
  # Volcano Plot
  p_vol <- ggplot2::ggplot(deg_results, ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=diff)) +
    ggplot2::geom_point(alpha=0.4) +
    ggplot2::scale_color_manual(values=c("Down"="blue", "Stable"="grey", "Up"="red")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 P-value")
  ggplot2::ggsave(file.path(output_dir, "1_Volcano_Plot.png"), p_vol, width = 8, height = 6, dpi = 300)
  # PCA Plot
  vsd <- DESeq2::vst(dds, blind=FALSE)
  p_pca <- DESeq2::plotPCA(vsd, intgroup = target_col) +
    ggplot2::theme_bw() + ggplot2::labs(title="PCA Plot")
  ggplot2::ggsave(file.path(output_dir, "2_PCA_Plot.png"), p_pca, width = 8, height = 6, dpi = 300)
  # Barplot Top 20
  top20 <- deg_results[order(deg_results$padj), ] |> head(20)
  p_bar <- ggplot2::ggplot(top20, ggplot2::aes(x=reorder(gene_symbol, log2FoldChange), y=log2FoldChange, fill=diff)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(title="Top 20 DEGs", x="Genes", y="Log2 Fold Change")
  ggplot2::ggsave(file.path(output_dir, "3_Top20_Barplot.png"), p_bar, width = 8, height = 8, dpi = 300)
  message(">>> Analysis completed successfully! Folder: ", output_dir)

  return(deg_results_output)
}
