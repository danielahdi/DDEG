#' DDEG: Advanced Differential Expression & Co-expression Analyzer (High-quality + Visualization)
#' @importFrom stats as.formula reorder
#' @importFrom utils head read.csv read.table write.csv
#' @importFrom graphics text abline par
#' @importFrom grDevices dev.off png recordPlot
#' @param count_file Path to the count matrix (CSV or TXT)
#' @param meta_file (Optional) Path to the metadata file
#' @param group_column (Optional) The column name in metadata for grouping
#' @param study_design (Optional) Vector of column names in metadata for module-trait analysis
#' @param p_val_cutoff Adjusted p-value cutoff (default 0.05)
#' @param logfc_cutoff Absolute log2FC cutoff (default 1)
#' @param run_enrichment Logical: Run GO/KEGG/Reactome enrichment? (default TRUE)
#' @param run_WGCNA Logical: Run WGCNA analysis? (default TRUE)
#' @param top_genes_wgcna Number of top variable genes for WGCNA (default 10000)
#' @param use_all_filtered_genes Logical: Use all filtered genes for WGCNA (default FALSE)
#' @param maxBlockSize Maximum genes per block for WGCNA (default 15000)
#' @export
run_DDEG <- function(count_file, meta_file = NULL, group_column = NULL,
                     study_design = NULL, p_val_cutoff = 0.05, logfc_cutoff = 1,
                     run_enrichment = TRUE, run_WGCNA = TRUE,
                     top_genes_wgcna = 10000, use_all_filtered_genes = FALSE,
                     maxBlockSize = 15000) {

  # -------------------------------
  # Load required packages
  # -------------------------------
  required_pkgs <- c("DESeq2","WGCNA","clusterProfiler","DOSE","ReactomePA",
                     "org.Hs.eg.db","AnnotationDbi","ggplot2","matrixStats","dplyr","SummarizedExperiment","parallel")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, update=FALSE, ask=FALSE)
  }

  library(DESeq2); library(WGCNA); library(clusterProfiler); library(DOSE)
  library(ReactomePA); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2); library(matrixStats); library(dplyr)
  library(SummarizedExperiment); library(parallel)

  options(stringsAsFactors = FALSE)
  enableWGCNAThreads(nThreads = parallel::detectCores()-1)

  # -------------------------------
  # Read count data
  # -------------------------------
  read_data <- function(path){
    if (grepl("\\.csv$", path)) read.csv(path, row.names=1, check.names=FALSE)
    else read.table(path, header=TRUE, row.names=1, check.names=FALSE)
  }
  counts <- read_data(count_file)
  counts <- counts[rowSums(counts) >= 10, ]
  rownames(counts) <- gsub("\\..*","",rownames(counts))
  counts <- round(counts)

  # -------------------------------
  # Read metadata
  # -------------------------------
  if (is.null(meta_file)){
    metadata <- data.frame(condition=factor(rep("Group1",ncol(counts))), row.names=colnames(counts))
    if (!is.null(group_column)) metadata[[group_column]] <- factor(metadata[[group_column]])
  } else {
    metadata <- read_data(meta_file)
    target_col <- if(is.null(group_column)) colnames(metadata)[1] else group_column
    common_samples <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[,common_samples]
    metadata <- metadata[common_samples, , drop=FALSE]
  }

  # -------------------------------
  # DESeq2 Differential Expression
  # -------------------------------
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=metadata,
                                        design = as.formula(paste0("~", ifelse(is.null(group_column), colnames(metadata)[1], group_column))))
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, independentFiltering=FALSE, cooksCutoff=FALSE)
  deg_results <- as.data.frame(res)

  # Annotation
  deg_results$gene_symbol <- mapIds(org.Hs.eg.db, keys=rownames(deg_results),
                                    column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  deg_results$gene_symbol[is.na(deg_results$gene_symbol)] <- rownames(deg_results)[is.na(deg_results$gene_symbol)]
  deg_results$log2FoldChange[is.na(deg_results$log2FoldChange)] <- 0
  deg_results$padj[is.na(deg_results$padj)] <- 1
  deg_results$diff <- "Stable"
  deg_results$diff[deg_results$log2FoldChange > logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Up"
  deg_results$diff[deg_results$log2FoldChange < -logfc_cutoff & deg_results$padj < p_val_cutoff] <- "Down"

  # -------------------------------
  # Output folder
  # -------------------------------
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- paste0("DDEG_Results_", timestamp)
  dir.create(output_dir, showWarnings=FALSE)
  deg_results$ensembl_id <- rownames(deg_results)
  write.csv(deg_results[, c("ensembl_id","gene_symbol","baseMean","log2FoldChange","pvalue","padj","diff")],
            file.path(output_dir,"DEG_Full_Table.csv"), row.names=FALSE)

  # VST for WGCNA
  vsd <- DESeq2::vst(dds, blind=FALSE)

  # -------------------------------
  # WGCNA Analysis
  # -------------------------------
  if (run_WGCNA){
    datExpr0 <- t(assay(vsd))
    message(">>> WGCNA: total filtered genes: ", ncol(datExpr0))

    if (use_all_filtered_genes) datExpr <- datExpr0
    else {
      rv <- colVars(datExpr0)
      select <- order(rv,decreasing=TRUE)[seq_len(min(top_genes_wgcna, ncol(datExpr0)))]
      datExpr <- datExpr0[,select]
      message(">>> Using top ", length(select), " variable genes")
    }

    # Remove outlier samples
    gsg <- goodSamplesGenes(datExpr, verbose=0)
    if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

    # Soft-thresholding
    sft <- pickSoftThreshold(datExpr, powerVector=1:20, networkType="signed", verbose=0)
    softPower <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)
    message(">>> Selected power: ", softPower)

    # Module detection
    net <- blockwiseModules(datExpr, power=softPower, TOMType="signed",
                            minModuleSize=20, mergeCutHeight=0.2,
                            maxBlockSize=maxBlockSize, verbose=0)
    colors <- labels2colors(net$colors)
    unique_colors <- unique(colors[colors!="grey"])
    message(">>> Detected ", length(unique_colors), " meaningful modules")

    # Module-trait correlation
    if (!is.null(study_design)){
      traitData <- metadata[,study_design, drop=FALSE]
      MEs <- net$MEs
      moduleTraitCor <- cor(MEs, traitData, use="p")
      moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
      write.csv(moduleTraitCor, file.path(output_dir,"WGCNA_ModuleTraitCor.csv"))
      write.csv(moduleTraitP, file.path(output_dir,"WGCNA_ModuleTraitPval.csv"))
    }

    # Hub genes
    MEs <- net$MEs
    moduleMembership <- bicor(datExpr, MEs, use="pairwise.complete.obs")
    hub_info <- data.frame()
    for (color in unique_colors){
      module_genes <- colnames(datExpr)[colors==color]
      mm <- abs(moduleMembership[module_genes,color])
      hub_info <- rbind(hub_info, data.frame(Gene=module_genes, Module=color, kME=mm))
    }
    top_hubs <- hub_info %>% group_by(Module) %>% slice_max(kME, n=10)
    write.csv(top_hubs, file.path(output_dir,"WGCNA_TopHubs_kME_PerModule.csv"), row.names=FALSE)

    saveRDS(net, file.path(output_dir,"WGCNA_Net.rds"))
    write.csv(data.frame(Gene=colnames(datExpr), ModuleColor=colors),
              file.path(output_dir,"WGCNA_Modules.csv"), row.names=FALSE)

    # -------------------------------
    # WGCNA Visualization
    # -------------------------------
    # Dendrogram + module colors
    png(file.path(output_dir,"WGCNA_Dendrogram.png"), width=1200, height=800)
    plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels=FALSE, hang=0.03,
                        addGuide=TRUE, guideHang=0.05)
    dev.off()

    # Module-trait heatmap
    if (!is.null(study_design)){
      textMatrix <- paste(signif(moduleTraitCor,2), "\n(", signif(moduleTraitP,1),")", sep="")
      png(file.path(output_dir,"WGCNA_ModuleTrait_Heatmap.png"), width=1200, height=1000)
      par(mar=c(6,8,3,3))
      labeledHeatmap(Matrix=moduleTraitCor,
                     xLabels=colnames(traitData),
                     yLabels=colnames(moduleTraitCor),
                     ySymbols=colnames(moduleTraitCor),
                     colorLabels=TRUE,
                     colors=blueWhiteRed(50),
                     textMatrix=textMatrix,
                     setStdMargins=FALSE,
                     cex.text=0.7,
                     zlim=c(-1,1),
                     main="Module-Trait Relationships")
      dev.off()
    }

    # Top hub gene barplots
    for (color in unique_colors){
      top_hubs_mod <- top_hubs %>% filter(Module==color)
      if(nrow(top_hubs_mod)>0){
        p <- ggplot(top_hubs_mod, aes(x=reorder(Gene,kME), y=kME)) +
          geom_bar(stat="identity", fill=color) + coord_flip() +
          theme_minimal() + labs(title=paste("Top hub genes - Module", color), y="kME", x="Gene")
        ggsave(file.path(output_dir,paste0("WGCNA_TopHubs_",color,".png")), p, width=8, height=6, dpi=300)
      }
    }
  }

  # -------------------------------
  # Enrichment Analysis
  # -------------------------------
  if (run_enrichment){
    deg_genes <- deg_results$ensembl_id[deg_results$diff!="Stable"]
    message(">>> Found ", length(deg_genes), " DEGs for enrichment")

    if(length(deg_genes)>0){
      # GO BP/MF/CC
      for (ont in c("BP","MF","CC")){
        ego <- enrichGO(deg_genes, OrgDb=org.Hs.eg.db, keyType="ENSEMBL",
                        ont=ont, pAdjustMethod="BH", qvalueCutoff=0.05)
        write.csv(as.data.frame(ego), file.path(output_dir,paste0("Enrichment_DEG_GO_",ont,".csv")), row.names=FALSE)
      }

      # KEGG & Reactome
      entrez <- bitr(deg_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
      if(length(entrez)>0){
        ekegg <- enrichKEGG(entrez, pAdjustMethod="BH", pvalueCutoff=0.05)
        ereact <- enrichPathway(entrez, pvalueCutoff=0.05)
        write.csv(as.data.frame(ekegg), file.path(output_dir,"Enrichment_DEG_KEGG.csv"), row.names=FALSE)
        write.csv(as.data.frame(ereact), file.path(output_dir,"Enrichment_DEG_Reactome.csv"), row.names=FALSE)
      }
    }

    # Module enrichment (GO BP)
    if (run_WGCNA && length(unique_colors)>0){
      for (color in unique_colors){
        module_genes <- colnames(datExpr)[colors==color]
        if(length(module_genes)<15) next
        ego_mod <- enrichGO(module_genes, OrgDb=org.Hs.eg.db, keyType="ENSEMBL",
                            ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05)
        if(nrow(as.data.frame(ego_mod))>0){
          write.csv(as.data.frame(ego_mod),
                    file.path(output_dir,paste0("Enrichment_Module_",color,"_GO_BP.csv")),
                    row.names=FALSE)

          # Dotplot
          png(file.path(output_dir,paste0("Module_",color,"_GO_BP_dotplot.png")), width=1200, height=800)
          print(dotplot(ego_mod, showCategory=15) + ggtitle(paste("GO BP - Module", color)))
          dev.off()
        }
      }
    }
  }

  message(">>> Analysis completed. Results in folder: ", output_dir)
  return(deg_results)
}

