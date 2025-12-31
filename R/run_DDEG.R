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
  message(">>> Pre-filtering: Removing genes with total counts < 10...")
  counts <- counts[rowSums(counts) >= 10, ]
  message(">>> Cleaning Gene IDs and preparing counts...")
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
  message(paste(">>> Running DESeq2 analysis grouped by:", target_col))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                        design = as.formula(paste0("~", target_col)))
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
  deg_results <- as.data.frame(res)
  message(">>> Mapping Ensembl IDs to Gene Symbols...")
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
  final_cols <- c("ensembl_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "diff")
  write.csv(deg_results[, final_cols], file.path(output_dir, "DEG_Full_Table.csv"), row.names = FALSE)
  # Volcano Plot
  p_vol <- ggplot2::ggplot(deg_results, ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=diff)) +
    ggplot2::geom_point(alpha=0.4) +
    ggplot2::scale_color_manual(values=c("Down"="blue", "Stable"="grey", "Up"="red")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-value")
  ggplot2::ggsave(file.path(output_dir, "1_Volcano_Plot.png"), p_vol, width = 8, height = 6, dpi = 300)
  # PCA Plot
  vsd <- DESeq2::vst(dds, blind=FALSE)
  p_pca <- DESeq2::plotPCA(vsd, intgroup = target_col) +
    ggplot2::theme_bw() + ggplot2::labs(title="PCA Plot")
  ggplot2::ggsave(file.path(output_dir, "2_PCA_Plot.png"), p_pca, width = 8, height = 6, dpi = 300)
  # Top 20 Barplot
  top20 <- deg_results[order(deg_results$padj), ] |> head(20)
  p_bar <- ggplot2::ggplot(top20, ggplot2::aes(x=reorder(gene_symbol, log2FoldChange), y=log2FoldChange, fill=diff)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(title="Top 20 DEGs", x="Genes", y="Log2 Fold Change")
  ggplot2::ggsave(file.path(output_dir, "3_Top20_Barplot.png"), p_bar, width = 8, height = 8, dpi = 300)
  # --- WGCNA Section ---
  cat("\nWould you like to run WGCNA analysis? (yes/no): ")
  user_response <- readline()
  if (tolower(trimws(user_response)) %in% c("yes", "y")) {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      BiocManager::install(c("WGCNA", "impute", "preprocessCore"), update = FALSE)
    }
    WGCNA::allowWGCNAThreads()
    cor <- WGCNA::cor  # رفع conflict
    datExpr <- t(SummarizedExperiment::assay(vsd))
    gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK) {
      datExpr <- datExpr[, gsg$goodGenes]
      message(">>> Removed problematic genes/samples based on goodSamplesGenes.")
    }
    # Sample clustering for outliers
    sampleTree <- hclust(dist(datExpr), method = "average")
    png(file.path(output_dir, "4_Sample_Clustering_Outliers.png"), width = 800, height = 600)
    plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    dev.off()
    message(">>> Check 4_Sample_Clustering_Outliers.png for potential sample outliers.")
    # Pick soft threshold
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 5, networkType = "signed")
    # SFT plot
    sft_df <- data.frame(Power = sft$fitIndices[,1],
                         SFT.R.sq = sft$fitIndices[,2],
                         Slope = sft$fitIndices[,3])
    p_sft <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power)) +
      ggplot2::geom_point(ggplot2::aes(y = SFT.R.sq), color = "red") +
      ggplot2::geom_text(ggplot2::aes(y = SFT.R.sq, label = Power), nudge_y = 0.05) +
      ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed", color = "blue") +
      ggplot2::labs(title = "Scale-Free Topology Fit", y = "Scale-Free R^2", x = "Soft Threshold Power") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(output_dir, "5_WGCNA_SFT_Plot.png"), p_sft, width = 8, height = 6)
    # Choose power
    softPower <- sft$powerEstimate
    if (is.na(softPower)) {
      softPower <- min(which(sft$fitIndices$SFT.R.sq >= 0.85))
      if (is.na(softPower)) softPower <- 9
    }
    message(">>> Selected soft threshold power: ", softPower)
    # Blockwise modules (signed network, large block size to minimize blocks)
    net <- WGCNA::blockwiseModules(datExpr, power = softPower,
                                   TOMType = "signed", networkType = "signed",
                                   minModuleSize = 30, verbose = 3, maxBlockSize = 25000)
    # Save network object and module assignments
    saveRDS(net, file.path(output_dir, "WGCNA_Network.rds"))
    gene_modules <- data.frame(Gene = colnames(datExpr),
                               Module = net$colors,
                               ModuleColor = WGCNA::labels2colors(net$colors))
    write.csv(gene_modules, file.path(output_dir, "WGCNA_Gene_Modules.csv"), row.names = FALSE)
    # Dendrogram plots per block
    n_blocks <- length(net$dendrograms)
    for (i in seq_len(n_blocks)) {
      block_genes <- net$blockGenes[[i]]
      colors_block <- net$colors[block_genes]
      png(file.path(output_dir, paste0("6_WGCNA_Dendrogram_Block", i, ".png")), width = 1000, height = 600)
      WGCNA::plotDendroAndColors(net$dendrograms[[i]],
                                 colors = WGCNA::labels2colors(colors_block),
                                 rowText = NULL,
                                 groupLabels = "Module Colors",
                                 dendroLabels = FALSE, hang = 0.03, addGuide = TRUE)
      dev.off()
    }
    # Module-trait correlations
    message(">>> Available metadata columns: ", paste(colnames(metadata), collapse = ", "))
    trait_col <- readline(prompt = "Enter column name for trait correlation: ")
    if (trimws(trait_col) != "" && trait_col %in% colnames(metadata)) {
      trait_data <- as.numeric(as.factor(metadata[[trait_col]]))
      names(trait_data) <- rownames(metadata)
      moduleTraitCor <- WGCNA::cor(net$MEs, trait_data, use = "p")
      moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))
      png(file.path(output_dir, "7_WGCNA_Module_Trait_Heatmap.png"), width = 600, height = 800)
      par(mar = c(6, 8.5, 3, 3))
      WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                            xLabels = trait_col,
                            yLabels = names(net$MEs),
                            ySymbols = names(net$MEs),
                            colorLabels = FALSE,
                            colors = WGCNA::blueWhiteRed(50),
                            textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                                               signif(moduleTraitP, 1), ")", sep = ""),
                            setStdMargins = FALSE,
                            cex.text = 0.8,
                            zlim = c(-1, 1),
                            main = paste("Module-Trait Relationships"))
      dev.off()
    }
  }
  message(">>> Analysis completed! Results saved in folder: ", output_dir)
  return(deg_results)
}
