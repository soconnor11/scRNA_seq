#docker run -it -v '/home/soconnor/old_home/ccNN/testData:/files' cplaisier/ccafv2_extra

#--------------------------------
# Set up section / load packages
#--------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library("org.Hs.eg.db")
library(aricode)
library(reticulate)
use_python('/usr/bin/python3')

# Set working directory
setwd("files/")

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)
mgenes = read.csv("ccAFv2_genes_full_dataset_102423_fc_25_v2_layers_600_200.csv")[,2]

# Extra
gene_ids = c('gene_symbols', 'ensembl')

# Gene sets
mesenchymal_geneset = c("CD44", "SPARCL1", "MKI67", "AURKA")
proneural_geneset = c("OLIG1", "OLIG2", "SOX4", "SOX8", "SOX2")
hypoxia_geneset = c("VEGFA", "ANGPTL4", "NDRG1", "CD44")
G0_G1_p53_targets = c("MDM2", "BBC3", "ZMAT3", "TP53I3", "CDKN1A")
goi_lst2 = list(mesenchymal_geneset, proneural_geneset, hypoxia_geneset, G0_G1_p53_targets)
sapply(goi_lst2, length)

neuralGO_geneset = c("SCRG1", "PLP1", "S100B", "GPM6B", "BEX1", "PTPRZ1", "PRCP", "PTN", "SOX4", "SAT1")
G1_geneset = c("IGFBP3", "IGFBP5", "MIAT", "MAP3K7CL", "AHNAK2", "TPST2", "DLG1" , "CMTM7", "C6orf15", "GJB2")
late_G1_geneset = c("EDN1", "CYR61", "ANKRD1", "CTGF", "PLK2", "UGCG", "ARID5B", "PLAU", "CCL2")
S_geneset = c("CCNE2", "CLSPN", "GINS2", "PCNA", "ATAD2", "MCM7", "MCM3", "SLBP", "GMNN", "KIAA0101")
s_g2_geneset = c("HIST1H4C", "CDK1", "HIST1H1E", "HIST1H1B", "UBE2C", "RRM2", "ZWINT", "HIST1H1C", "HMGB2")
G2_M_geneset = c("CCNB1", "CENPF", "CKS2", "PTTG1", "CDC20", "TOP2A", "NUSAP1", "CENPA")
M_early_G1_geneset = c("HMGN2", "TUBA1B", "STMN1", "BIRC5", "HMGB1", "TROAP", "HNRNPA2B1", "H2AFZ", "ARL6IP1")
goi_lst = list(neuralGO_geneset,G1_geneset, late_G1_geneset, S_geneset, s_g2_geneset, G2_M_geneset, M_early_G1_geneset )
sapply(goi_lst, length)

# SCDataHugnotIGF features
hugnot = c('S100B', 'SOX2', 'SOX4', 'MKI67', 'APOE', 'VIM', 'CLU', 'FABP7','OLIG1','OLIG2', 'DLL3', 'HES6')
# convert to ensembl IDs
ensembl_hugnot = mapIds(org.Hs.eg.db, keys = hugnot, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_hugnot = na.omit(data.frame(ensembl_hugnot))
ensembl_hugnot_plot = ensembl_hugnot$ensembl_hugnot

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# convert to ensembl IDs
ensembl_s_genes = mapIds(org.Hs.eg.db, keys = s.genes, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_s_genes_2 = na.omit(data.frame(ensembl_s_genes))$ensembl_s_genes
ensembl_g2m_genes = mapIds(org.Hs.eg.db, keys = g2m.genes, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_g2m_genes_2 = na.omit(data.frame(ensembl_g2m_genes))$ensembl_g2m_genes


# Plotting order & colors
ccAF_colors = c("G1/other" = "#9aca3c", "Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
ccAF_order = c("G1/other", 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccSeurat_order = c('G1', 'S', 'G2M')


roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}


# Plot multibar heatmap
suppressPackageStartupMessages({
  library(rlang)
})
DoMultiBarHeatmap <- function (object,
                               features = NULL,
                               cells = NULL,
                               group.by = "ident",
                               additional.group.by = NULL,
                               additional.group.sort.by = NULL,
                               cols.use = NULL,
                               group.bar = TRUE,
                               disp.min = -2.5,
                               disp.max = NULL,
                               slot = "scale.data",
                               assay = NULL,
                               label = TRUE,
                               size = 5.5,
                               hjust = 0,
                               angle = 45,
                               raster = TRUE,
                               draw.lines = TRUE,
                               lines.width = NULL,
                               group.bar.height = 0.02,
                               combine = TRUE)
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }

  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ",
                paste(bad.sorts, collapse = ", "))
      }
    }
  }

  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                             slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }

    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]

    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) *
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }

      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }

      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }

    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])

    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster,
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features,
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])

    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))

        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }

        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])

        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(plot +
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off"))

        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos),
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity",
                                   data = label.x.pos, aes_string(label = "group",
                                                                  x = "label.x.pos"), y = y.max + y.max *
                                     0.03 * 0.5, angle = angle, hjust = hjust,
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) *
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}


# Process data individually
scProcessData = function(res_dir, tag, assay = 'SCT', cutoff = 0.5, species= 'human', resolution = 0.8, save_dir = 'analysis_output', obj_dir = 'seurat_objects', analysis_dir = 'cutoff_analysis', symbol = F){
  cat('\n',tag,'\n')
  # Set up folders
  dir.create(file.path(res_dir, 'FINAL'), showWarnings = FALSE)
  resdir1 = file.path(res_dir, tag)
  resdir2 = file.path(resdir1, save_dir)
  resdir3 = file.path(resdir1, obj_dir)
  resdir4 = file.path(res_dir, 'FINAL')
  #---------------------------
  # Load filtered / normalized data
  #---------------------------
  gene_id = 'ensembl'
  seurat2 = readRDS(file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id),'.rds')))
  #seurat2 = readRDS(file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id),'.rds')))
  # Load ccSeurat calls
  ccseurat_calls = read.csv(file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')), row.names = 'X')
  seurat2 <- AddMetaData(seurat2, ccseurat_calls, col.name = "Phase")
  sub4 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub4])
  # Load ccAF calls
  ccAF_calls = read.csv(file.path(resdir2, paste0(tag, '_ccAF_calls.csv')), row.names = 'X')
  seurat2 <- AddMetaData(seurat2, ccAF_calls$ccAF, col.name = "ccAF")
  sub2 = ccAF_order %in% factor(seurat2$ccAF)
  seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAF_order[sub2])
  #---------------------------
  # Classify with ccAFv2
  #---------------------------
  seurat2 = PredictCellCycle(seurat2, assay = assay, cutoff = cutoff, species=species, gene_id=gene_id)
  sub3 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub3])
  tmp = data.frame(table(seurat2$ccAFv2))
  rownames(tmp) = tmp$Var1
  print((tmp['Freq']/dim(seurat2)[2])*100)
  write.csv(seurat2$ccAFv2, file.path(resdir4, paste0(tag, '_ccAFv2_calls.csv')))
  write.csv(data.frame((tmp['Freq']/dim(seurat2)[2])*100), file.path(resdir4, paste0(tag, '_ccAFv2_call_frequency.csv')))
  #---------------------------
  # Normalize
  #---------------------------
  seurat2 = SCTransform(seurat2, verbose = FALSE)
  # Add marker gene counts for each cell
  seurat_subset = seurat2[mgenes]
  mgene_counts = GetAssayData(seurat_subset, slot = "counts")
  non_zero_mgenes = colSums(mgene_counts > 0)
  cat('Cells that have non-zero ccAFv2 genes: ', length(non_zero_mgenes), '\n')
  # Add as meta data column
  seurat2 <- AddMetaData(seurat2, non_zero_mgenes, col.name = 'ccAFv2_mgene_counts')
  seurat2 <- RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 <- RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  Idents(seurat2) <- seurat2$ccAFv2
  # Find cluster marker genes
  cluster_markers= FindAllMarkers(seurat2, only.pos = TRUE)
  cluster_markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 0.5) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
  cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
  cluster_markers$gene = cluster_markers_genes
  write.csv(cluster_markers, file.path(resdir4, paste0(tag,"_scTransform_Markers_together.csv")))
  top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('Plotting UMAPs and ccAFv2 barplot \n')
  d1 = DimPlot(seurat2, reduction = "umap", label=F, group.by="seurat_clusters", raster = FALSE) + ggtitle("seurat_clusters")
  d2 = DimPlot(seurat2, reduction = "umap", label=F, group.by="Phase", cols = ccSeurat_colors, raster = FALSE)  + ggtitle("Phase")
  d3 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAF", cols = ccAF_colors[sub2]) + ggtitle("ccAF")
  d4 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAFv2", cols = ccAFv2_colors[sub3]) + ggtitle("ccAFv2")
  pdf(file.path(resdir4, paste0(tag, '.pdf')), width = 10, height = 8)
  lst = list(d1, d2, d3, d4)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = "")
  # gene sets
  p1 = lapply(ensembl_hugnot_plot, function(goi) {
    if(goi %in% rownames(seurat2)){
      fp1 = FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_hugnot)[ensembl_hugnot$ensembl_hugnot == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
      return(fp1 + FontSize(x.title = 10, y.title = 10))
    }
  })
  grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3, 4), c(5,6,7,8), c(9,10,11, 12)), top = "")
  #------- ccAFv2 vs. seurat clusters - cell percentages -------#
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(table(seurat2$ccAFv2, seurat2$seurat_clusters), beside = FALSE, col = ccAFv2_colors[sub3], xlab = "clusters ID", ylab = "Cell count", legend.text = rownames(table(seurat2$ccAFv2, seurat2$seurat_clusters)), args.legend=list(title="ccAFv2 classification"))
  #--- ccAFv2 vs. cluster ids stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$seurat_clusters)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$seurat_clusters)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$seurat_clusters)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 4.5, y=max(colSums(cf_1)), bty = "n"))
  #--- ccAFv2 vs. ccSeurat stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$Phase)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$Phase)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$Phase)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 1.5, y=max(colSums(cf_1)), bty = "n"))
  #--- ccAFv2 and number of marker genes boxplot ---#
  v1 = VlnPlot(seurat2, features = 'ccAFv2_mgene_counts', group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + theme(legend.position = 'none') + xlab("ccAFv2") + ylab("ccAFv2 marker gene counts")
  v2 = VlnPlot(seurat2, features = 'nCount_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + theme(legend.position = 'none') + xlab("ccAFv2") + ylab("nCount_RNA")
  v3 = VlnPlot(seurat2, features = 'nFeature_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub1]) + theme(legend.position = 'none') + xlab("ccAFv2") + ylab("nFeature_RNA")
  lst2 = list(v2, v3, v1)
  grid.arrange(grobs = lst2, layout_matrix = rbind(c(1, 2), c(3, NA)), top = "")
  print(DoHeatmap(object = seurat2, features = names(top10_symbol), group.colors = ccAFv2_colors[sub1], size = 4) + scale_y_discrete(labels = top10_symbol))
  dev.off()
  # Save processed data out as loom object
  # Change factored metadata to characters
  seurat2$ccAF = as.character(seurat2$ccAF)
  seurat2$ccAFv2 = as.character(seurat2$ccAFv2)
  seurat2$Phase = as.character(seurat2$Phase)
  cat('saving processed data as loom and rds...\n')
  data_loom_2 <- as.loom(seurat2, file.path(resdir4, paste0(tag, '_processed.loom')), verbose = FALSE, overwrite = TRUE)
  data_loom_2$close_all()
  saveRDS(seurat2, file.path(resdir4, paste0(tag, '_processed.rds')))
  return(seurat2)
}


plotComparison = function(res_dir, tag, conditions = c('GF', 'noGF'), save_dir = 'analysis_output'){
  resdirs = res_dir
  tags = tag
  extratags = conditions
  for(resdir in resdirs){
    new_tags = c()
    for(tag1 in tags){
      for(extratag1 in extratags){
        new_tags = append(new_tags, paste0(tag1, '_', extratag1))
        }
      for(new_tag in new_tags){
        cat('\n',new_tag,'\n')
        #resdir1 = file.path(resdir, tag, paste0(new_tag))
        resdir1 = file.path(resdir, tag, 'FINAL')
        #resdir2 = file.path(resdir1, save_dir)
        data = read.csv(file.path(resdir1, paste0(new_tag, '_ccAFv2_call_frequency.csv')), row.names = 'X')
        # Set up matrices
        if(new_tag == new_tags[1]){
          gf_data = matrix(nrow = length(rownames(data)), ncol = 3)
          gf_data[,1] = data$Freq
          gf_data[,2] = new_tag
          gf_data[,3] = rownames(data)
        } else if(new_tag == new_tags[2]){
          nogf_data = matrix(nrow = length(rownames(data)), ncol = 3)
          nogf_data[,1] = data$Freq
          nogf_data[,2] = new_tag
          nogf_data[,3] = rownames(data)
          }
        }
      }
    }
    tmp = rbind(gf_data, nogf_data)
    df = data.frame(tmp)
    df$X1 = as.numeric(df$X1)
    colnames(df) <- c('percentage','condition', 'ccAFv2_state')
    df$percentage = round(df$percentage, 2)
    pdf(file.path(resdir, tag, 'FINAL/ccAFv2_comparison_between_conditions.pdf'), width = 10, height = 8)
    print(ggplot(df, aes(x = condition, y = percentage, fill = factor(ccAFv2_state, levels = ccAFv2_order), label = paste0(percentage,'%'))) + geom_bar(position='stack', stat='identity') + geom_text(
      size = 3, position = position_stack(vjust = 0.5)) + scale_fill_manual(values=ccAFv2_colors, name='ccAFv2') + theme(plot.margin = margin(1,1,1.5,1.2, "cm")))
    dev.off()
  }


# Integrate and process data
scIntegrateAndProcessData = function(res_dir, tag, conditions = c('GF', 'noGF'), resolution = 0.6, save_dir = 'integration_output', obj_dir = 'seurat_objects', analysis_dir = 'analysis_output', symbol = F){
  cat('\n',tag,'\n')
  # Set up folders
  resdir1 = file.path(res_dir, tag)
  dir.create(file.path(resdir1, 'FINAL'), showWarnings = FALSE)
  resdir4 = file.path(resdir1, 'FINAL')
  #dir.create(file.path(resdir1, save_dir))
  #resdir2 = file.path(resdir1, save_dir)
  #dir.create(file.path(resdir2, obj_dir))
  #resdir3 = file.path(resdir2, obj_dir)
  #dir.create(file.path(resdir2, analysis_dir))
  #---------------------------
  # Load normalized data
  #---------------------------
  gene_id = 'ensembl'
  datas = list()
  cat('Load data\n')
  for(cond1 in conditions){
    datas[[cond1]] = readRDS(file.path(resdir1, paste0(tag, '_', cond1), obj_dir, paste0(tag, '_', cond1, '_filtered_', gene_id,'.rds')))
    # Specify condition in object
    datas[[cond1]]$condition = cond1
    datas[[cond1]]
    # Load ccSeurat calls
    ccseurat_calls = read.csv(file.path(resdir1, paste0(tag, '_',cond1), analysis_dir, paste0(tag, '_', cond1, '_ccSeurat_calls.csv')), row.names = 'X')
    datas[[cond1]] = AddMetaData(datas[[cond1]], ccseurat_calls, col.name = "Phase")
    # Load ccAF calls
    ccAF_calls = read.csv(file.path(resdir1, paste0(tag, '_',cond1), analysis_dir, paste0(tag, '_',cond1, '_ccAF_calls.csv')), row.names = 'X')
    datas[[cond1]] = AddMetaData(datas[[cond1]], ccAF_calls$ccAF, col.name = "ccAF")
    # Load ccAFv2 calls
    ccAFv2_calls = read.csv(file.path(resdir4, paste0(tag, '_',cond1, '_ccAFv2_calls.csv')), row.names = 'X')
    datas[[cond1]] = AddMetaData(datas[[cond1]], ccAFv2_calls$x, col.name = "ccAFv2")
  }
  # merge and integrate
  cat('Merge and integrate GF and noGF datasets\n')
  seurat1 <- merge(datas[["GF"]], datas[["noGF"]], add.cell.ids = c("GF", "noGF"))
  seurat.list <- SplitObject(seurat1, split.by = "condition")
  seurat.list <- lapply(X = seurat.list, FUN = SCTransform, verbose= FALSE)
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",anchor.features = features)
  seurat.combined.sct <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")
  seurat2 = seurat.combined.sct
  # Order ccSeurat calls
  sub4 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub4])
  # Order ccAF calls
  sub2 = ccAF_order %in% factor(seurat2$ccAF)
  seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAF_order[sub2])
  # Order ccAFv2 calls
  sub3 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub3])
  tmp = data.frame(table(seurat2$ccAFv2, seurat2$condition))
  write.csv(tmp, file.path(resdir4, paste0(tag, '_integrated_ccAFv2_counts_per_conditions.csv')))
  write.csv(seurat2$ccAFv2, file.path(resdir4, paste0(tag, '_integrated_ccAFv2_calls.csv')))
  seurat2 <- RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 <- RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  Idents(seurat2) <- seurat2$ccAFv2
  # Find cluster marker genes
  cluster_markers= FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE)
  cluster_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
  cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
  cluster_markers$gene = cluster_markers_genes
  write.csv(cluster_markers, file.path(resdir4, paste0(tag, "_integrated_scTransform_Markers_together.csv")))
  top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('Plotting UMAPs and ccAFv2 barplots \n')
  d1 = DimPlot(seurat2, reduction = "umap", label=F, group.by="seurat_clusters", split.by = 'condition') + ggtitle("seurat_clusters")
  d2 = DimPlot(seurat2, reduction = "umap", label=F, group.by="Phase", cols = ccSeurat_colors, split.by = 'condition')  + ggtitle("Phase")
  d3 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAF", cols = ccAF_colors[sub2], split.by = 'condition') + ggtitle("ccAF")
  d4 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAFv2", cols = ccAFv2_colors[sub3], split.by = 'condition') + ggtitle("ccAFv2")
  # ccAFv2 barplot
  tmp2 = tmp[tmp$Var2 == 'GF',]
  tmp3 = tmp2 %>% mutate(perc=Freq/sum(Freq))
  tmp4 = tmp[tmp$Var2 == 'noGF',]
  tmp5 = tmp4 %>% mutate(perc=Freq/sum(Freq))
  df = rbind(tmp3, tmp5)
  # cluster barplot
  tmp6 = data.frame(table(seurat2$seurat_clusters, seurat2$ccAFv2, seurat2$condition))
  tmp7 = tmp6[tmp6$Var3=='GF',]
  yend7 = roundUpNice(max(tmp7$Freq))
  tmp8 = tmp6[tmp6$Var3=='noGF',]
  yend8 = roundUpNice(max(tmp8$Freq))
  max1 = max(yend7, yend8)
  pdf(file.path(resdir4, paste0(tag, '_integrated.pdf')), width = 12, height = 8)
  lst = list(d1, d2, d3, d4)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, 4)), top = "")
  print(ggplot(data = df, aes(x=Var1, y=perc, fill=Var2)) + geom_bar(stat="identity", position = "dodge") + scale_y_continuous(labels = scales::percent) + scale_fill_discrete(name = "Condition") + xlab("ccAFv2") + ylab("Frequency"))
  b1 = ggplot(data=tmp7, aes(x=Var1, y=Freq,  fill=Var2)) + geom_bar(position='stack', stat='identity') + xlab("seurat_clusters") + ylab("cells") + ggtitle('GF') + scale_fill_manual(values=ccAFv2_colors, name='ccAFv2') + ylim(0, max1)
  b2 = ggplot(data=tmp8, aes(x=Var1, y=Freq,  fill=Var2)) + geom_bar(position='stack', stat='identity') + xlab("seurat_clusters") + ylab("cells") + ggtitle('noGF') + scale_fill_manual(values=ccAFv2_colors, name='ccAFv2') + ylim(0, max1)
  lst2 = list(b1, b2)
  grid.arrange(grobs = lst2, layout_matrix = cbind(c(1, 2)), top = "")
  #print(DoHeatmap(object = seurat2, features = names(top10_symbol), group.colors = ccAFv2_colors[sub3], size = 3, angle = 45) + scale_y_discrete(labels = top10_symbol))
  cols.use <- list(ccAFv2=ccAFv2_colors)
  print(DoMultiBarHeatmap(object = seurat2, features = names(top10_symbol), group.by = 'ccAFv2', additional.group.by= "condition", size = 3, angle = 45, cols.use=cols.use) + scale_y_discrete(labels = top10_symbol) + scale_color_manual(values=ccAFv2_colors[sub3]))
  dev.off()
  # Save processed data out as loom object
  # Change factored metadata to characters
  seurat2$ccAF = as.character(seurat2$ccAF)
  seurat2$ccAFv2 = as.character(seurat2$ccAFv2)
  seurat2$Phase = as.character(seurat2$Phase)
  cat('saving processed data as loom and rds...\n')
  data_loom_2 <- as.loom(seurat2, file.path(resdir4, paste0(tag, '_integrated_processed.loom')), verbose = FALSE, overwrite = TRUE)
  data_loom_2$close_all()
  saveRDS(seurat2, file.path(resdir4, paste0(tag, '_integrated_processed.rds')))
  return(seurat2)
}


#---------------------------------------------------
# Downstream analysis & plotting (all together)
#--------------------------------------------------

# SCDataHugnotIGF
hu_ensembl_process = list()
hu_ensembl_process[['LGG275_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_GF', resolution = 0.61)
hu_ensembl_process[['LGG275_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_noGF', resolution = 0.15)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'LGG275', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'LGG275', conditions = c('GF', 'noGF'), resolution = 0.3)

hu_ensembl_process[['LGG85_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_GF', resolution = 0.6)
hu_ensembl_process[['LGG85_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'LGG85', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'LGG85', conditions = c('GF', 'noGF'), resolution = 0.4)

hu_ensembl_process[['LGG336_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_GF', resolution = 0.6)
hu_ensembl_process[['LGG336_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_noGF', resolution = 0.4)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'LGG336', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'LGG336', conditions = c('GF', 'noGF'), resolution = 0.3)

hu_ensembl_process[['LGG349_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_GF', resolution = 0.6)
hu_ensembl_process[['LGG349_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'LGG349', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'LGG349', conditions = c('GF', 'noGF'), resolution = 0.2)

hu_ensembl_process[['BT054_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_GF', resolution = 0.6)
hu_ensembl_process[['BT054_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'BT054', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'BT054', conditions = c('GF', 'noGF'), resolution = 0.4)

hu_ensembl_process[['BT088_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_GF', resolution = 0.6)
hu_ensembl_process[['BT088_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'BT088', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'BT088', conditions = c('GF', 'noGF'), resolution = 0.2)

hu_ensembl_process[['BT138_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_GF', resolution = 0.6)
hu_ensembl_process[['BT138_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'BT138', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'BT138', conditions = c('GF', 'noGF'), resolution = 0.3)

hu_ensembl_process[['BT237_GF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_GF', resolution = 0.6)
hu_ensembl_process[['BT237_noGF']] = scProcessData(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_noGF', resolution = 0.6)
plotComparison(res_dir = 'SCDataHugnotIGF', tag = 'BT237', conditions = c('GF', 'noGF'))
scIntegrateAndProcessData(res_dir = 'SCDataHugnotIGF', tag = 'BT237', conditions = c('GF', 'noGF'), resolution = 0.3)
