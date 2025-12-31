#' DDEG: Differential Data Expression Generator with WGCNA
#' @importFrom stats as.formula reorder
#' @importFrom utils head read.csv read.table write.csv
#' @importFrom graphics text abline par
#' @importFrom grDevices dev.off png recordPlot
#' @param count_file Path to the count matrix (CSV or TXT)
#' @param meta_file (Optional) Path to the metadata file
#' @param group_column (Optional) The column name in metadata
#' @param p_val_cutoff Cutoff for adjusted p-value (default 0.05)
#' @param logfc_cutoff Cutoff for absolute log2 fold change (default 1)
#' @export
run_DDEG <- function(count_file, meta_file = NULL, group_column = NULL, p_val_cutoff = 0.05, logfc_cutoff = 1) {
  # رفع Note مربوط به متغیرهای گلوبال
  log2FoldChange <- padj <- gene_symbol <- diff <- NULL
  message(">>> Loading libraries and data...")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  }
  read_data <- function(path) {
    if (grepl("\\.csv", path)) return(read.csv(path, row.names = 1, check.names = FALSE))
    return(read.table(path, header = TRUE, row.names = 1, check.names = FALSE))
  }
  counts <- read_data(count_file)
  counts <- counts[rowSums(counts) >= 10, ]
  rownames(counts) <- gsub("\\..*", "", rownames(counts))
  counts <- round(counts)
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
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                        design = as.formula(paste0("~", target_col)))
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
  deg_results <- as.data.frame(res)
  deg_results$gene_symbol <- tryCatch({
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(deg_results),
                          column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  }, error = function(e) { rep(NA, nrow(deg_results)) })
  deg_results$gene_symbol[is.na(deg_results$gene_symbol)] <- rownames(deg_results)[is.na(deg_results$gene_symbol)]
  deg_results$log2FoldChange[is.na(deg_results$log2FoldChange)] <- 0
  deg_results$pvalue[is.na(deg_results$pvalue)] <- 1
  deg_results$padj[is.na(deg_results$padj)] <- 1
  deg_results$diff <- "Stable"
  deg_results$diff[deg_results$log2FoldChange > logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Up"
  deg_results$diff[deg_results$log2FoldChange < -logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Down"
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- paste0("DDEG_Results_", timestamp)
  if (!dir.exists(output_dir)) dir.create(output_dir)
  deg_results$ensembl_id <- rownames(deg_results)
  write.csv(deg_results[, c("ensembl_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "diff")],
            file.path(output_dir, "DEG_Full_Table.csv"), row.names = FALSE)
  # Plots
  p_vol <- ggplot2::ggplot(deg_results, ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=diff)) +
    ggplot2::geom_point(alpha=0.4) + ggplot2::scale_color_manual(values=c("Down"="blue", "Stable"="grey", "Up"="red")) +
    ggplot2::theme_minimal() + ggplot2::labs(title="Volcano Plot")
  ggplot2::ggsave(file.path(output_dir, "1_Volcano_Plot.png"), p_vol)
  vsd <- DESeq2::vst(dds, blind=FALSE)
  p_pca <- DESeq2::plotPCA(vsd, intgroup = target_col) + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(output_dir, "2_PCA_Plot.png"), p_pca)
  # --- WGCNA Section ---
  cat("\nWould you like to run WGCNA analysis? (yes/no): ")
  user_response <- readline()
  if (tolower(trimws(user_response)) %in% c("yes", "y")) {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      BiocManager::install(c("WGCNA", "impute", "preprocessCore"), update = FALSE)
    }
    WGCNA::allowWGCNAThreads()
    # جدید: رفع کانفلیکت cor
    cor <- WGCNA::cor
    datExpr <- t(SummarizedExperiment::assay(vsd))
    gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
    if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes]
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = c(1:20), verbose = 0)
    softPower <- ifelse(is.na(sft$powerEstimate), 9, sft$powerEstimate)
    net <- WGCNA::blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                                   minModuleSize = 30, verbose = 0)
    png(file.path(output_dir, "4_WGCNA_Dendrogram.png"), width = 800, height = 600)
    WGCNA::plotDendroAndColors(net$dendrograms[[1]], WGCNA::labels2colors(net$colors),
                               "Module Colors", dendroLabels = FALSE)
    dev.off()
    message(">>> Available metadata columns: ", paste(colnames(metadata), collapse=", "))
    trait_col <- readline(prompt = "Enter column name for correlation: ")
    if (trait_col %in% colnames(metadata)) {
      trait_data <- as.numeric(as.factor(metadata[[trait_col]]))
      moduleTraitCor <- WGCNA::cor(net$MEs, trait_data, use = "p")
      moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))
      png(file.path(output_dir, "5_WGCNA_Heatmap.png"), width = 600, height = 800)
      par(mar = c(6, 8, 3, 3))
      WGCNA::labeledHeatmap(Matrix = moduleTraitCor, xLabels = trait_col,
                            yLabels = names(net$MEs), colors = WGCNA::blueWhiteRed(50),
                            textMatrix = paste0(round(moduleTraitCor, 2), "\n(", round(moduleTraitPvalue, 3), ")"),
                            cex.text = 0.8, main = "Module-Trait Relationships")
      dev.off()
    }
  }
  return(deg_results)
}
