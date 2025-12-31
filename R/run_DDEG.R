#' DDEG: Differential Data Expression Generator with WGCNA (Optimized for Speed)
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
  # --- Optimized WGCNA Section ---
  cat("\nWould you like to run WGCNA analysis? (yes/no): ")
  user_response <- readline()
  if (tolower(trimws(user_response)) %in% c("yes", "y")) {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      BiocManager::install(c("WGCNA", "impute", "preprocessCore"), update = FALSE)
    }
    # فعال کردن multithreading
    WGCNA::enableWGCNAThreads(nThreads = parallel::detectCores() - 1)  # استفاده حداکثری از هسته‌ها
    cor <- WGCNA::cor

    datExpr0 <- t(SummarizedExperiment::assay(vsd))

    # بهینه‌سازی اصلی: فیلتر به top 5000-8000 ژن متغیر (این کار زمان TOM رو از ساعت‌ها به دقیقه کم می‌کنه)
    nTop <- 5000  # می‌تونی به 5000 کم کنی اگر هنوز کند بود
    rv <- matrixStats::rowVars(datExpr0)
    select <- order(rv, decreasing = TRUE)[seq_len(min(nTop, length(rv)))]
    datExpr <- datExpr0[, select]
    message(paste0(">>> Filtered to top ", nTop, " most variable genes for fast WGCNA analysis."))

    gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

    # Sample outlier detection
    sampleTree <- hclust(dist(datExpr), method = "average")
    png(file.path(output_dir, "4_Sample_Clustering.png"), width = 1000, height = 600)
    plot(sampleTree, main = "Sample Clustering (Outliers)", xlab = "", sub = "")
    dev.off()

    # Pick soft threshold (signed برای حفظ همبستگی منفی)
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 6:20, networkType = "signed", verbose = 5)

    # SFT Plot
    p_sft <- ggplot2::ggplot(data.frame(Power = sft$fitIndices[,1], R2 = sft$fitIndices[,2]),
                             ggplot2::aes(x = Power, y = R2)) +
      ggplot2::geom_point(color = "red") + ggplot2::geom_text(ggplot2::aes(label = Power), nudge_y = 0.05) +
      ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed") +
      ggplot2::theme_minimal() + ggplot2::labs(title = "Scale-Free Topology Fit")
    ggplot2::ggsave(file.path(output_dir, "5_WGCNA_SFT_Plot.png"), p_sft)

    softPower <- sft$powerEstimate
    if (is.na(softPower)) softPower <- 12  # fallback رایج برای signed

    message(">>> Using power: ", softPower)

    # BlockwiseModules با تنظیمات سریع
    net <- WGCNA::blockwiseModules(datExpr, power = softPower,
                                   TOMType = "signed", networkType = "signed",
                                   minModuleSize = 30, mergeCutHeight = 0.25,
                                   verbose = 3, maxBlockSize = 15000,  # بلوک کوچک برای سرعت
                                   nThreads = parallel::detectCores() - 1)

    saveRDS(net, file.path(output_dir, "WGCNA_Network.rds"))
    gene_modules <- data.frame(Gene = colnames(datExpr),
                               ModuleColor = WGCNA::labels2colors(net$colors))
    write.csv(gene_modules, file.path(output_dir, "WGCNA_Gene_Modules.csv"), row.names = FALSE)

    # Dendrogram per block
    for (i in seq_len(length(net$dendrograms))) {
      block_genes <- net$blockGenes[[i]]
      png(file.path(output_dir, paste0("6_Dendrogram_Block", i, ".png")), width = 1200, height = 700)
      WGCNA::plotDendroAndColors(net$dendrograms[[i]],
                                 WGCNA::labels2colors(net$colors[block_genes]),
                                 groupLabels = "Module Colors", dendroLabels = FALSE)
      dev.off()
    }

    # Trait correlation
    message(">>> Metadata columns: ", paste(colnames(metadata), collapse = ", "))
    trait_col <- readline(prompt = "Enter trait column (or Enter to skip): ")
    if (trimws(trait_col) != "" && trait_col %in% colnames(metadata)) {
      trait_data <- as.numeric(as.factor(metadata[[trait_col]]))
      moduleTraitCor <- WGCNA::cor(net$MEs, trait_data, use = "p")
      moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))
      png(file.path(output_dir, "7_Module_Trait_Heatmap.png"), width = 700, height = 900)
      WGCNA::labeledHeatmap(Matrix = moduleTraitCor, xLabels = trait_col, yLabels = names(net$MEs),
                            colors = WGCNA::blueWhiteRed(50),
                            textMatrix = paste(signif(moduleTraitCor, 2), "\n(p=", signif(moduleTraitP, 1), ")", sep=""),
                            cex.text = 0.9, main = "Module-Trait Relationships")
      dev.off()
    }
  }
  message(">>> Analysis completed! Folder: ", output_dir)
  return(deg_results)
}
