# Run docker
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
library(ccAFv2)
library("org.Hs.eg.db")
library(aricode)

# Set working directory
setwd("files/")

# Extra
gene_ids = c('gene_symbols', 'ensembl')

# Mitochondrial genes as ensembl IDs
mito_genes = read_csv("mito_genes.csv", show_col_types = FALSE) %>% pull(mito)

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


#------------------------------------------------------
# Functions
#---------------------------------------------------

scQC = function(res_dir, tag, mito_genes, v1 = 2000, v2 = 30000, h1 = 0.001, h2 = 0.1, data_dir = 'outs', save_dir = 'analysis_output', obj_dir = 'seurat_objects', analysis_dir = 'cutoff_analysis',  mt = 'MT-', symbol = F){
  cat('\n',tag,'\n')
  resdir1 = file.path(res_dir, tag)
  # Create folders
  dir.create(file.path(resdir1, save_dir), showWarnings = FALSE)
  resdir2 = file.path(resdir1, save_dir)
  #dir.create(file.path(save_dir), showWarnings = FALSE)
  #resdir2 = file.path(save_dir)
  dir.create(file.path(resdir1, obj_dir), showWarnings = FALSE)
  resdir3 = file.path(resdir1, obj_dir)
  dir.create(file.path(resdir1, analysis_dir), showWarnings = FALSE)
  resdir4 = file.path(resdir1, analysis_dir)
  #---------------------
  # Load in data
  #---------------------
  gene_column = 1
  gene_id = 'ensembl'
  if(symbol){
    gene_column = 2
    gene_id = 'gene_symbols'
  }
  data = Read10X(file.path(resdir1, data_dir), gene.column=gene_column) # column 1 is ensembl (in 10X mtx file)
  cat('Raw data', dim(data)[2], 'cells', dim(data)[1], 'genes \n')
  # Substitute underscores if necessary
  rownames(data) = gsub("_", "-", rownames(data))
  # Create seurat object
  seurat1 = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
  cat('Basic filter', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
  if(symbol){
    mito_genes = grep(mt, rownames(seurat1))
  }
  seurat1[['percent.mito']] = PercentageFeatureSet(seurat1, features = mito_genes)/100
  #---------------------
  # Quality control
  #---------------------
  cat('Quality control \n')
  # Quality control plots for choosing cutoffs
  pdf(file.path(resdir2, paste0(tag, '_QC_plot_to_choose_cutoffs.pdf')))
  plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  # Quality control filtering
  keep.detect = which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 & seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)
  seurat1 = subset(seurat1, cells=colnames(seurat1)[keep.detect])
  # Save out filtered object
  cat('Filtered to', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
  saveRDS(seurat1, file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id), '.rds')))
  # Use for ccAF application in python
  data_loom <- as.loom(seurat1, file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id), '.loom')), verbose = FALSE, overwrite = TRUE)
  data_loom$close_all()
  # Save as new object so can go back to previous non-normalized / scaled seurat object if need too
  seurat2 = seurat1
  #------------------------------------------------------
  # Normalization with sctransform
  #---------------------------------------------------
  cat('Normalization \n')
  seurat2 = SCTransform(seurat2, verbose = FALSE)
  cat('Normalized genes:', dim(seurat2@assays$SCT@data)[1], 'features,', length(seurat2@assays$SCT@var.features), 'highly variable genes \n')
  # Classify with ccSeurat and save out as csv
  if(symbol){
    seurat2 <- CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
    write.csv(seurat2$Phase, file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')))
  }
  cat('Saving normalized RDS object \n')
  saveRDS(seurat2, file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id), '.rds')))
  return(seurat2)
}


#--------------------------------------
# Quality control & data preparation
#--------------------------------------

# SCDataHugnotIGF
hu_ensembl = list()
hu_ensembl[['LGG275_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_GF', mito_genes = mito_genes, v1 = 5000, v2 = 76000, h1 = 0.001, h2 = 0.15)
hu_ensembl[['LGG275_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 35000, h1 = 0.001, h2 = 0.15)

hu_ensembl[['LGG336_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_GF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.01, h2 = 0.15)
hu_ensembl[['LGG336_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 140000, h1 = 0.01, h2 = 0.15)

hu_ensembl[['LGG85_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_GF', mito_genes = mito_genes, v1 = 6000, v2 = 80000, h1 = 0.01, h2 = 0.15)
hu_ensembl[['LGG85_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.03, h2 = 0.17)

hu_ensembl[['LGG349_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_GF', mito_genes = mito_genes, v1 = 6000, v2 = 100000, h1 = 0.01, h2 = 0.10)
hu_ensembl[['LGG349_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.001, h2 = 0.08)

hu_ensembl[['BT054_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_GF', mito_genes = mito_genes, v1 = 6000, v2 = 150000, h1 = 0.02, h2 = 0.20)
hu_ensembl[['BT054_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.24)

hu_ensembl[['BT088_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_GF', mito_genes = mito_genes, v1 = 6000, v2 = 130000, h1 = 0.01, h2 = 0.12)
hu_ensembl[['BT088_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 105000, h1 = 0.001, h2 = 0.10)

hu_ensembl[['BT138_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_GF', mito_genes = mito_genes, v1 = 6000, v2 = 75000, h1 = 0.01, h2 = 0.10)
hu_ensembl[['BT138_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_noGF', mito_genes = mito_genes, v1 = 6500, v2 = 93000, h1 = 0.01, h2 = 0.13)

hu_ensembl[['BT237_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_GF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.01, h2 = 0.15)
hu_ensembl[['BT237_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.17)


hu_symbol = list()
hu_symbol[['LGG275_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_GF', mito_genes = mito_genes, v1 = 5000, v2 = 76000, h1 = 0.001, h2 = 0.15, symbol = T)
hu_symbol[['LGG275_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG275', tag = 'LGG275_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 35000, h1 = 0.001, h2 = 0.15, symbol = T)

hu_symbol[['LGG336_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_GF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.01, h2 = 0.15, symbol = T)
hu_symbol[['LGG336_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG336', tag = 'LGG336_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 140000, h1 = 0.01, h2 = 0.15, symbol = T)

hu_symbol[['LGG85_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_GF', mito_genes = mito_genes, v1 = 6000, v2 = 80000, h1 = 0.01, h2 = 0.15, symbol = T)
hu_symbol[['LGG85_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG85', tag = 'LGG85_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.03, h2 = 0.17, symbol = T)

hu_symbol[['LGG349_GF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_GF', mito_genes = mito_genes, v1 = 6000, v2 = 100000, h1 = 0.01, h2 = 0.10, symbol = T)
hu_symbol[['LGG349_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/LGG349', tag = 'LGG349_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 120000, h1 = 0.001, h2 = 0.08, symbol = T)

hu_symbol[['BT054_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_GF', mito_genes = mito_genes, v1 = 6000, v2 = 150000, h1 = 0.02, h2 = 0.20, symbol = T)
hu_symbol[['BT054_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT054', tag = 'BT054_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.24, symbol = T)

hu_symbol[['BT088_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_GF', mito_genes = mito_genes, v1 = 6000, v2 = 130000, h1 = 0.01, h2 = 0.12, symbol = T)
hu_symbol[['BT088_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT088', tag = 'BT088_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 105000, h1 = 0.001, h2 = 0.10, symbol = T)

hu_symbol[['BT138_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_GF', mito_genes = mito_genes, v1 = 6000, v2 = 75000, h1 = 0.01, h2 = 0.10, symbol = T)
hu_symbol[['BT138_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT138', tag = 'BT138_noGF', mito_genes = mito_genes, v1 = 6500, v2 = 93000, h1 = 0.01, h2 = 0.13, symbol = T)

hu_symbol[['BT237_GF']] = scQC(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_GF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.01, h2 = 0.15, symbol = T)
hu_symbol[['BT237_noGF']] = scQC(res_dir = 'SCDataHugnotIGF/BT237', tag = 'BT237_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.17, symbol = T)
