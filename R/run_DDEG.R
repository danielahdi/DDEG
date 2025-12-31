#' DDEG: Advanced Differential Expression & Co-expression Analyzer (Balanced Speed & Quality)
#' @importFrom stats as.formula reorder
#' @importFrom utils head read.csv read.table write.csv
#' @importFrom graphics text abline par
#' @importFrom grDevices dev.off png recordPlot
#' @param count_file Path to the count matrix (CSV or TXT)
#' @param meta_file (Optional) Path to the metadata file
#' @param group_column (Optional) The column name in metadata for grouping
#' @param p_val_cutoff Adjusted p-value cutoff (default 0.05)
#' @param logfc_cutoff Absolute log2FC cutoff (default 1)
#' @param run_enrichment Logical: Run GO/KEGG/Reactome enrichment? (default TRUE)
#' @param top_genes_wgcna Number of top variable genes for WGCNA (default 15000 for better quality)
#' @param use_all_filtered_genes Logical: Use ALL filtered genes (no variance filter) for maximum quality? (default FALSE - slower)
#' @export
run_DDEG <- function(count_file, meta_file = NULL, group_column = NULL,
                     p_val_cutoff = 0.05, logfc_cutoff = 1,
                     run_enrichment = TRUE, top_genes_wgcna = 15000,
                     use_all_filtered_genes = FALSE) {
  # متغیرهای گلوبال
  log2FoldChange <- padj <- gene_symbol <- diff <- NULL

  message(">>> Starting HIGH-QUALITY DDEG analysis (balanced speed & biological meaning)...")

  # لود کتابخانه‌ها
  required_pkgs <- c("org.Hs.eg.db", "clusterProfiler", "DOSE", "enrichplot", "ggplot2", "matrixStats", "dplyr", "WGCNA", "parallel", "ReactomePA")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    library(pkg, character.only = TRUE)
  }

  # خواندن داده
  read_data <- function(path) {
    if (grepl("\\.csv$", path)) read.csv(path, row.names = 1, check.names = FALSE)
    else read.table(path, header = TRUE, row.names = 1, check.names = FALSE)
  }
  counts <- read_data(count_file)
  message(">>> Initial genes: ", nrow(counts))
  counts <- counts[rowSums(counts) >= 10, ]
  message(">>> After low-count filter: ", nrow(counts), " genes")
  rownames(counts) <- gsub("\\..*", "", rownames(counts))
  counts <- round(counts)

  # metadata
  if (is.null(meta_file)) {
    col_names <- colnames(counts)
    group_labels <- rep("Treatment", length(col_names))
    control_idx <- grep("control|ctrl|normal|healthy|wild|wt", col_names, ignore.case = TRUE)
    if (length(control_idx) == 0) {
      half <- ceiling(length(col_names)/2)
      group_labels[1:half] <- "Group1"; group_labels[(half+1):length(col_names)] <- "Group2"
    } else group_labels[control_idx] <- "Control"
    metadata <- data.frame(condition = factor(group_labels), row.names = col_names)
    target_col <- "condition"
  } else {
    metadata <- read_data(meta_file)
    target_col <- if (is.null(group_column)) colnames(metadata)[1] else group_column
    common_samples <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, common_samples]
    metadata <- metadata[common_samples, ]
  }

  # DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                        design = as.formula(paste0("~", target_col)))
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
  deg_results <- as.data.frame(res)

  # Annotation
  deg_results$gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                   keys = rownames(deg_results),
                                                   column = "SYMBOL", keytype = "ENSEMBL",
                                                   multiVals = "first")
  deg_results$gene_symbol[is.na(deg_results$gene_symbol)] <- rownames(deg_results)[is.na(deg_results$gene_symbol)]

  # NA handling & diff status
  deg_results$log2FoldChange[is.na(deg_results$log2FoldChange)] <- 0
  deg_results$pvalue[is.na(deg_results$pvalue)] <- 1
  deg_results$padj[is.na(deg_results$padj)] <- 1
  deg_results$diff <- "Stable"
  deg_results$diff[deg_results$log2FoldChange > logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Up"
  deg_results$diff[deg_results$log2FoldChange < -logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Down"

  # خروجی و پلات‌های پایه
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- paste0("DDEG_HighQuality_Results_", timestamp)
  dir.create(output_dir, showWarnings = FALSE)
  deg_results$ensembl_id <- rownames(deg_results)
  write.csv(deg_results[, c("ensembl_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "diff")],
            file.path(output_dir, "DEG_Full_Table.csv"), row.names = FALSE)

  # Plots
  p_vol <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(padj), color = diff)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Down" = "blue", "Stable" = "grey", "Up" = "red")) +
    theme_minimal() + labs(title = "Volcano Plot")
  ggsave(file.path(output_dir, "1_Volcano.png"), p_vol, width = 8, height = 6, dpi = 300)

  vsd <- DESeq2::vst(dds, blind = FALSE)
  p_pca <- DESeq2::plotPCA(vsd, intgroup = target_col) + theme_bw()
  ggsave(file.path(output_dir, "2_PCA.png"), p_pca, width = 8, height = 6, dpi = 300)

  top20 <- head(deg_results[order(deg_results$padj), ], 20)
  p_bar <- ggplot(top20, aes(x = reorder(gene_symbol, log2FoldChange), y = log2FoldChange, fill = diff)) +
    geom_bar(stat = "identity") + coord_flip() + theme_light() + labs(title = "Top 20 DEGs")
  ggsave(file.path(output_dir, "3_Top20.png"), p_bar, width = 8, height = 8, dpi = 300)

  # WGCNA با تمرکز روی کیفیت (بیشتر ژن + تنظیمات برای moduleهای معنی‌دار)
  cat("\nRun WGCNA (high quality mode)? (yes/no): ")
  if (tolower(trimws(readline())) %in% c("yes", "y")) {
    enableWGCNAThreads(nThreads = parallel::detectCores() - 1)
    cor <- WGCNA::cor

    datExpr0 <- t(SummarizedExperiment::assay(vsd))
    message(">>> Total filtered genes for WGCNA: ", nrow(datExpr0))

    if (use_all_filtered_genes) {
      datExpr <- datExpr0
      message(">>> Using ALL filtered genes (maximum quality - may be slower)")
    } else {
      rv <- rowVars(datExpr0)
      select <- order(rv, decreasing = TRUE)[seq_len(min(top_genes_wgcna, nrow(datExpr0)))]
      datExpr <- datExpr0[, select]
      message(paste(">>> Using top", length(select), "most variable genes (balanced quality/speed)"))
    }

    gsg <- goodSamplesGenes(datExpr, verbose = 0)
    if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

    # power گسترده‌تر برای پیدا کردن بهترین fit
    sft <- pickSoftThreshold(datExpr, powerVector = 1:30, networkType = "signed", verbose = 0)
    softPower <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)
    message(">>> Selected power: ", softPower, " (for scale-free topology)")

    # تنظیمات برای moduleهای معنی‌دارتر
    net <- blockwiseModules(datExpr, power = softPower, TOMType = "signed",
                            minModuleSize = 20,  # کوچکتر برای moduleهای بیشتر
                            mergeCutHeight = 0.15,  # merge قوی‌تر برای moduleهای coherent
                            verbose = 0, maxBlockSize = 20000)

    colors <- net$colors
    unique_colors <- unique(colors[colors != "grey"])
    message(paste(">>> Detected", length(unique_colors), "meaningful modules"))

    if (length(unique_colors) == 0) {
      message(">>> No meaningful modules - check data variability or try use_all_filtered_genes = TRUE")
    } else {
      # Hub genes با kME (module membership) - دقیق‌تر از kWithin
      MEs <- net$MEs
      moduleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
      names(moduleMembership) <- gsub("ME", "", names(moduleMembership))
      hub_info <- data.frame()
      for (color in unique_colors) {
        module_genes <- colnames(datExpr)[colors == color]
        mm <- abs(moduleMembership[module_genes, color])
        hub_info <- rbind(hub_info, data.frame(Gene = module_genes, Module = color, kME = mm))
      }
      top_hubs <- hub_info %>% group_by(Module) %>% top_n(10, wt = kME)  # top 10 hubs
      write.csv(top_hubs, file.path(output_dir, "WGCNA_TopHubs_kME_PerModule.csv"), row.names = FALSE)
    }

    saveRDS(net, file.path(output_dir, "WGCNA_Net.rds"))
    write.csv(data.frame(Gene = colnames(datExpr), ModuleColor = labels2colors(colors)),
              file.path(output_dir, "WGCNA_Modules.csv"), row.names = FALSE)
  }

  # Enrichment کامل‌تر و معنی‌دار (GO BP/MF/CC + KEGG + Reactome)
  if (run_enrichment) {
    deg_genes <- deg_results$ensembl_id[deg_results$diff != "Stable"]
    message(paste(">>> Found", length(deg_genes), "DEGs for enrichment"))

    if (length(deg_genes) > 0) {
      ego_bp <- enrichGO(deg_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05)
      ego_mf <- enrichGO(deg_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pvalueCutoff = 0.05)
      ego_cc <- enrichGO(deg_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC", pvalueCutoff = 0.05)

      write.csv(as.data.frame(ego_bp), file.path(output_dir, "Enrichment_DEG_GO_BP.csv"), row.names = FALSE)
      write.csv(as.data.frame(ego_mf), file.path(output_dir, "Enrichment_DEG_GO_MF.csv"), row.names = FALSE)
      write.csv(as.data.frame(ego_cc), file.path(output_dir, "Enrichment_DEG_GO_CC.csv"), row.names = FALSE)

      if (nrow(as.data.frame(ego_bp)) > 0) {
        dotplot(ego_bp) + ggtitle("GO BP - DEGs")
        ggsave(file.path(output_dir, "Enrichment_DEG_GO_BP.png"))
      }

      # KEGG & Reactome
      entrez <- bitr(deg_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
      if (length(entrez) > 0) {
        ekegg <- enrichKEGG(entrez, pvalueCutoff = 0.05)
        ereact <- enrichPathway(entrez, pvalueCutoff = 0.05)
        write.csv(as.data.frame(ekegg), file.path(output_dir, "Enrichment_DEG_KEGG.csv"), row.names = FALSE)
        write.csv(as.data.frame(ereact), file.path(output_dir, "Enrichment_DEG_Reactome.csv"), row.names = FALSE)
      }
    }

    # Enrichment moduleها (فقط BP برای معنی‌داری بیولوژیکی)
    if (exists("net") && length(unique_colors) > 0) {
      for (color in unique_colors) {
        module_genes <- colnames(datExpr)[colors == color]
        if (length(module_genes) < 15) next
        ego_mod <- enrichGO(module_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05)
        if (nrow(as.data.frame(ego_mod)) > 0) {
          write.csv(as.data.frame(ego_mod), file.path(output_dir, paste0("Enrichment_Module_", color, "_GO_BP.csv")), row.names = FALSE)
          dotplot(ego_mod, showCategory = 15) + ggtitle(paste("GO BP - Module", color))
          ggsave(file.path(output_dir, paste0("Enrichment_Module_", color, "_GO_BP.png")))
        }
      }
    }
  }

  message(">>> High-quality analysis completed! More genes, better modules, meaningful hubs (kME), full enrichment.")
  message(">>> Tip: If few modules, try use_all_filtered_genes = TRUE for maximum coverage.")
  return(deg_results)
}
