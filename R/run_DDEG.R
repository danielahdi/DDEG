#' DDEG: Advanced Differential Expression & Co-expression Analyzer (Highly Optimized & Unique Features)
#' @importFrom stats as.formula reorder
#' @importFrom utils head read.csv read.table write.csv
#' @importFrom graphics text abline par
#' @importFrom grDevices dev.off png recordPlot
#' @param count_file Path to the count matrix (CSV or TXT)
#' @param meta_file (Optional) Path to the metadata file
#' @param group_column (Optional) The column name in metadata for grouping
#' @param p_val_cutoff Adjusted p-value cutoff (default 0.05)
#' @param logfc_cutoff Absolute log2FC cutoff (default 1)
#' @param run_enrichment Logical: Run GO/KEGG enrichment on DEGs and WGCNA modules? (default TRUE)
#' @param top_genes_wgcna Number of top variable genes for WGCNA (default 8000)
#' @export
run_DDEG <- function(count_file, meta_file = NULL, group_column = NULL,
                     p_val_cutoff = 0.05, logfc_cutoff = 1,
                     run_enrichment = TRUE, top_genes_wgcna = 8000) {
  # متغیرهای گلوبال
  log2FoldChange <- padj <- gene_symbol <- diff <- NULL

  message(">>> Starting advanced DDEG analysis...")

  # نصب/لود کتابخانه‌های ضروری
  required_pkgs <- c("org.Hs.eg.db", "clusterProfiler", "DOSE", "enrichplot", "ggplot2", "matrixStats")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }

  # خواندن داده
  read_data <- function(path) {
    if (grepl("\\.csv$", path)) read.csv(path, row.names = 1, check.names = FALSE)
    else read.table(path, header = TRUE, row.names = 1, check.names = FALSE)
  }
  counts <- read_data(count_file)
  counts <- counts[rowSums(counts) >= 10, ]
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
  output_dir <- paste0("DDEG_Advanced_Results_", timestamp)
  dir.create(output_dir, showWarnings = FALSE)
  deg_results$ensembl_id <- rownames(deg_results)
  write.csv(deg_results[, c("ensembl_id", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "diff")],
            file.path(output_dir, "DEG_Full_Table.csv"), row.names = FALSE)

  # Volcano, PCA, Top20
  library(ggplot2)
  ggplot(deg_results, aes(x=log2FoldChange, y=-log10(padj), color=diff)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("Down"="blue", "Stable"="grey", "Up"="red")) +
    theme_minimal() + labs(title="Volcano Plot")
  ggsave(file.path(output_dir, "1_Volcano.png"), width=8, height=6, dpi=300)

  vsd <- DESeq2::vst(dds, blind=FALSE)
  DESeq2::plotPCA(vsd, intgroup=target_col) + theme_bw()
  ggsave(file.path(output_dir, "2_PCA.png"), width=8, height=6, dpi=300)

  top20 <- head(deg_results[order(deg_results$padj), ], 20)
  ggplot(top20, aes(x=reorder(gene_symbol, log2FoldChange), y=log2FoldChange, fill=diff)) +
    geom_bar(stat="identity") + coord_flip() + theme_light()
  ggsave(file.path(output_dir, "3_Top20.png"), width=8, height=8, dpi=300)

  # WGCNA بهینه‌شده
  cat("\nRun WGCNA? (yes/no): ")
  if (tolower(trimws(readline())) %in% c("yes", "y")) {
    WGCNA::enableWGCNAThreads(nThreads = parallel::detectCores() - 1)
    cor <- WGCNA::cor

    datExpr0 <- t(SummarizedExperiment::assay(vsd))
    rv <- matrixStats::rowVars(datExpr0)
    select <- order(rv, decreasing = TRUE)[seq_len(min(top_genes_wgcna, nrow(datExpr0)))]
    datExpr <- datExpr0[, select]
    message(paste(">>> Using top", top_genes_wgcna, "variable genes for WGCNA"))

    gsg <- WGCNA::goodSamplesGenes(datExpr)
    if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 6:20, networkType = "signed", verbose = 0)
    softPower <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)

    net <- WGCNA::blockwiseModules(datExpr, power = softPower, TOMType = "signed",
                                   minModuleSize = 30, mergeCutHeight = 0.25,
                                   verbose = 0, maxBlockSize = 15000)

    saveRDS(net, file.path(output_dir, "WGCNA_Net.rds"))
    write.csv(data.frame(Gene = colnames(datExpr), Module = WGCNA::labels2colors(net$colors)),
              file.path(output_dir, "WGCNA_Modules.csv"), row.names = FALSE)

    # Hub genes per module
    hub_info <- data.frame()
    for (color in unique(net$colors)) {
      if (color == "grey") next
      module_genes <- colnames(datExpr)[net$colors == color]
      connectivity <- intramodularConnectivity(WGCNA::adjacency(datExpr[, module_genes]), net$colors[net$colors == color])
      hub_info <- rbind(hub_info, data.frame(Gene = module_genes, Module = color, kWithin = connectivity$kWithin))
    }
    top_hubs <- hub_info |> dplyr::group_by(Module) |> dplyr::top_n(5, kWithin)
    write.csv(top_hubs, file.path(output_dir, "WGCNA_TopHubs_PerModule.csv"), row.names = FALSE)
  }

  # Enrichment جدید و خاص (GO + KEGG + Reactome)
  if (run_enrichment) {
    library(clusterProfiler)
    deg_genes <- deg_results$ensembl_id[deg_results$diff != "Stable"]
    ego <- enrichGO(deg_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL")
    ekegg <- enrichKEGG(bitr(deg_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
    ereact <- enrichPathway(bitr(deg_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)

    write.csv(as.data.frame(ego), file.path(output_dir, "Enrichment_DEG_GO.csv"), row.names = FALSE)
    write.csv(as.data.frame(ekegg), file.path(output_dir, "Enrichment_DEG_KEGG.csv"), row.names = FALSE)
    write.csv(as.data.frame(ereact), file.path(output_dir, "Enrichment_DEG_Reactome.csv"), row.names = FALSE)

    # پلات‌های enrichment
    dotplot(ego) + ggtitle("GO Enrichment - DEGs")
    ggsave(file.path(output_dir, "Enrichment_GO_Dotplot.png"))

    # اگر WGCNA ران شد، enrichment روی moduleها
    if (exists("net")) {
      for (color in unique(net$colors[net$colors != "grey"])) {
        module_genes <- colnames(datExpr)[net$colors == color]
        ego_mod <- enrichGO(module_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
        if (!is.null(ego_mod)) {
          dotplot(ego_mod) + ggtitle(paste("GO BP - Module", color))
          ggsave(file.path(output_dir, paste0("Enrichment_Module_", color, "_GO.png")))
        }
      }
    }
  }

  message(">>> Advanced analysis completed! Unique features added: Hub genes + Multi-pathway enrichment.")
  return(deg_results)
}
