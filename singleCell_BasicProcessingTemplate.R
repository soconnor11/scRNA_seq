
#--------------------------------
# Set up section / load packages
#--------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
sessionInfo()
library(ssgsea.GBM.classification)
library(verification)
library(MCMCpack)


#------------------------------------------------------
# Read in data section / set up seurat object / QC
#---------------------------------------------------

# Set working directory
setwd("files/")
resdir = "redo_analysis"


#------------------------------------------------------
# Load in data
#---------------------------------------------------

tag = 'hBMMSC' # useful for saving figures
data_dir <- 'hBMMSC/filtered_feature_bc_matrix/' # where your data is located

#------ USE ENSEMBL IDs IF WANT TO SAVE OUT AS LOOM FILE TO APPLY CCAF --------#
# Genes as ensembl IDs - export as loom file to apply ccAF classification
#data <- Read10X(data.dir = data_dir, gene.column=1) # column 1 is ensembl IDs (in 10X mtx file)
#dim(data) #[1] 36601  7207
#------------------------------------------------------------------------------#

# Load in data
# Genes as gene symbols
data <- Read10X(data.dir = data_dir, gene.column=2) # column 2 is gene symbols (in 10X mtx file)
dim(data)

# Create seurat object
seurat1 <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
seurat1
dim(seurat1)


#------------------------------------------------------
# Quality control
#---------------------------------------------------

#------ USE ENSEMBL IDs IF WANT TO SAVE OUT AS LOOM FILE TO APPLY CCAF --------#
#mitogenes <- read_csv("hBMMSC/mito_genes.csv")
#list = as.list(mitogenes)
#seurat1[["percent.mito"]] <- PercentageFeatureSet(seurat1, features = as.list(mitogenes)$mito)/100
#------------------------------------------------------------------------------#

# Genes as gene symbols
# Find mitochondiral genes
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")

# Plot & adjust cutoffs before filtering seurat object
pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 6000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()

# Before filtering data - look at above plot and adjust cutoffs!
# Once find optimal cutoffs, enter values into keep.detect

# Quality control filtering
keep.detect <- which(seurat1@meta.data$percent.mito < 0.1 & seurat1@meta.data$percent.mito > 0.009 &
                        seurat1@meta.data$nCount_RNA < 130000 & seurat1@meta.data$nCount_RNA > 6000)
length(keep.detect) # will give you amount of cells after filtering
ncol(seurat1) - length(keep.detect) # will tell you how many cells you filtered out with above cutoffs
seurat1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect]) # actually filter
dim(seurat1) # gives you dimensions of filtered object

# Check if filtering worked - cells should be inside red line cutoffs (if outside, you have your < > backwards!)
# Make sure to match cutoff values with filtered values
pdf(file.path(resdir, paste0(tag, "_QC_plot_POST.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 6000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()


#------------------------------------------------------
# Save out as loom file (check to make sure genes are ensembl IDs)
# Comment out when finished
#---------------------------------------------------

#data_loom <- as.loom(seurat1, filename = paste0(tag, "_data.loom"), verbose = FALSE)
#data_loom
#data_loom$close_all()


#------------------------------------------------------
# Downstream processing
#---------------------------------------------------

# Load ccAF classifications and add as column into data
#ccAF <- read_csv("hBMMSC_data_ccAF_calls.csv")
#seurat1 <- AddMetaData(object = seurat1, metadata = ccAF$ccAF, col.name = 'ccAF')


#------------------------------------------------------
# Clustering / seurat pipeline with fastMNN -scTransform version.
#---------------------------------------------------

# save as new object so can go back to previous non-normalized / scaled seurat object if need too
seurat2 <- seurat1
#res = 0.8 # Choose resolution

seurat2 <- SCTransform(seurat2, vars.to.regress = "nCount_RNA", verbose = FALSE)
#seurat2 <- FindVariableFeatures(seurat2, selection.method = "vst", nfeatures = 2000)
seurat2 <- RunPCA(seurat2, dims = 1:30)
#seurat2 <- RunTSNE(seurat2, dims = 1:30)
seurat2 <- RunUMAP(seurat2,   dims = 1:30)
seurat2 <- FindNeighbors(seurat2,   dims = 1:30)
seurat2 <- FindClusters(seurat2)
#seurat2 <- FindClusters(seurat2, resolution = res)

pdf(file.path(resdir, paste0(tag, "_umap_test.pdf")), width = 8, height = 8)
m1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6) # if you use TNSE, replace "umap" with "tsne" in reduction parameter
dev.off()

# Look to see if number of clusters makes sense with your data. If not, un-comment line 128 and change resolution, and run up to 141. Looking at marker genes to find meaning for clusters also helps!
# Run through code until you find cluster_markers variable - literature search and then return in line 129 if want to change # of clusters

# Visualize
m1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6) + ggtitle ("de novo clustering")
#af = DimPlot(seurat2, reduction = "umap", group.by="ccAF", label=T) + ggtitle("ccAF")


#------------------------------------------------------
# Make Violin Plots for each cluster
#---------------------------------------------------

cluster_no = length(unique(seurat2$seurat_clusters))
o1 = VlnPlot(object = seurat2, features= "nFeature_RNA",   pt.size=0.1 )
o2 = VlnPlot(object = seurat2, features= "nCount_RNA",  pt.size=0.1 )
o4 = VlnPlot(object = seurat2, features= "percent.mito",   pt.size=0.1)


#---------------------------------------------------
# cell cycle scoring - ccSeurat
#---------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat2 <- CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

c1 =  DimPlot(seurat2, reduction = "umap", label=T, group.by="Phase") + ggtitle("ccSeurat")

#---------------------------------------------------
# Find marker genes and save to excel file.
#---------------------------------------------------

log2_thres = log2(1.25)
fdr_thres=0.05
cluster_markers= FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE) # positive marker gens
table(cluster_markers$cluster)
cluster_no = unique(cluster_markers$cluster)

sp = split(cluster_markers, cluster_markers$cluster)
names(sp) = paste0("cluster-", names(sp))
write_xlsx(sp, path = file.path(resdir, paste0(tag, "_scTransform_Markers_by_cluster.xlsx")))

cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10


#---------------------------------------------------
# save as RDS object and write plots to pdf
#---------------------------------------------------
saveRDS(seurat2 , file.path(resdir, paste0(tag, "_seurat2.rds")))

# Dump out pdf of plots
pdf(paste0(tag, "_test.pdf"), width =8, height = 10)

# umaps ( with cell cycle scoring )
#lst = list(m1, c1, af )
#grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,NA), c(3,NA)), top = "")

lst = list(m1, c1 ) # no ccAF
grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,NA)), top = "")

lst2 = list( o1, o2, o4)
grid.arrange(grobs = lst2, layout_matrix = rbind(c(1,2), c(3)), top = "")

# Heatmap
DoHeatmap(object = seurat2, features = top10$gene)

dev.off()


#---------------------------------------------------
# Assigning cell type identity to clusters
#---------------------------------------------------

new.cluster.ids <- c("celltype1", "celltype2", "celltype3") # number of new cluster IDs needs to be same as number of leiden clusters (also in same order - leiden cluster 0 = "celltype1" in this list and so on...)
names(new.cluster.ids) <- levels(seurat2)
seurat2 <- RenameIdents(seurat2, new.cluster.ids) # rename clusters in seurat object
m2 = DimPlot(seurat2, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + ggtitle("labeled clusters")

# Plot labeled umap
pdf(file.path(resdir, paste0(tag, "_labeled_umap.pdf")), width =8, height = 8)
m2
dev.off()


#---------------------------------------------------
# Additional resources!!!
#---------------------------------------------------

# use Seurat API and PBMC tutorial for help with paramater functions & plot examples!
# https://satijalab.org/seurat/
