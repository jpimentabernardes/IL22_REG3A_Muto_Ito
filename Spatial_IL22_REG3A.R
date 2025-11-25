---
  title: "Spatial RNAseq Go"
output: html_document
date: "2025-09-04"
---
  
library(Seurat)     # Main package for single-cell and spatial transcriptomics analysis
library(SeuratObject)
library(patchwork)  # Helps combine multiple ggplot2 plots into one figure
library(tidyverse)  # Essential collection of R packages, including ggplot2, dplyr, magrittr
library(grid)    # Plot multiple objects
library(viridis) # A series of color maps that are designed to improve graph readability
library(SeuratWrappers) #SeuratWrappers is a collection of community-provided methods and extensions for Seurat.
library(Matrix)
library(readr)
library(dplyr)
library(tidyr)
library(arrow)  # for parquet files
library(jsonlite)
library(OpenImageR)


#library(spacexr) #learning cell types and cell type-specific differential expression in spatial transcriptomics data.
# For other plots

# Colors
cal_pal50 = c("#Fa1a8e", "#009B7D", "#ff9933", "#7083a4", "#ffce45", "#015e05", 
              "#fedf25", "#d2b48c", "#bb55e1", "#6ec5ab", "#5d57af", "#143341", 
              "#761445", "#d65b5a", "#94043a", "#e7a6cd", "#204519", "#87afd1", 
              "#9b9a9d", "#f95b70", "#83c874", "#808080", "#452b5b", "#ecb100", 
              "#f46124", "#525252", "#4c84a3", "#00bfff", "#01b4c6", "#174d79", 
              "#a6a0f2", "#76facc", "#8491be", "#a32a2f", "#1c8859", "#2cc012", 
              "#35782b", "#9c6755", "#3b3960", "#eeb0a1", "#3e1e13", "#0064c3", 
              "#d81e4a", "#74646c", "#f675da", "#ffce45", "#ec7014", "#e50000", 
              "#000000", "#a4527c", "#041859")

setwd('/Users/joana/Documents/IKMB/spatialRNAseq/Go_IL22_spatial/')

B21 = Load10X_Spatial("input/B21_07755/outs/binned_outputs/square_008um/", assay="RNA")
B21$orig.ident = "B21"
Idents(B21) = "B21"
B21

regionB21<-read.csv('input/B21_07755/outs/binned_outputs/square_008um/region.csv')

# Make sure the barcodes are rownames in regionB21
rownames(regionB21) <- regionB21$Barcode
# Subset regionB21 to only include barcodes present in B21
regionB21 <- regionB21[colnames(B21), , drop = FALSE]
# Add the region column to B21 metadata
B21$region <- regionB21$region
# Quick check
head(B21@meta.data)
table(B21$region)


B22 = Load10X_Spatial("input/B22_07388/outs/binned_outputs/square_008um/", assay="RNA")
B22$orig.ident = "B22"
Idents(B22) = "B22"
B22

regionB22<-read.csv('input/B22_07388/outs/binned_outputs/square_008um/region.csv')

rownames(regionB22) <- regionB22$Barcode
regionB22 <- regionB22[colnames(B22), , drop = FALSE]
B22$region <- regionB22$region
head(B22@meta.data)
table(B22$region)

# Add sample identifiers (important before merging)
Idents(B21) <- "B21"
Idents(B22) <- "B22"
B21 <- RenameCells(B21, add.cell.id = "B21")
B22 <- RenameCells(B22, add.cell.id = "B22")

# Visualize to confirm
SpatialDimPlot(
  B21,
  group.by = "region",
  cols = c("SI" = "skyblue", "Colon" = "tomato"),
  pt.size.factor = 1.6
) 

SpatialDimPlot(
  B22,
  group.by = "region",
  cols = c("SI" = "skyblue", "Colon" = "tomato"),
  pt.size.factor = 1.6
) 

# Merge
combined <- merge(B21, y = B22)

# Check result
combined
table(combined$orig.ident)

Idents(combined) <- "orig.ident"


## 2. Quality control

p1 = VlnPlot(B21, group.by = 'orig.ident',features="nCount_RNA", pt.size=0, layer="counts") + NoLegend()
p2 = VlnPlot(B22, group.by = 'orig.ident',features="nCount_RNA", pt.size=0, layer="counts") + NoLegend()
p1 | p2


p1 = SpatialFeaturePlot(B21, features="nCount_RNA")
p2 = SpatialFeaturePlot(B22, features="nCount_RNA")
p1 | p2

p1 = VlnPlot(combined, group.by = 'orig.ident',features="nCount_RNA", pt.size=0, layer="counts") + NoLegend()
p2 = SpatialFeaturePlot(combined, features="nCount_RNA")
p1 | p2

p1 = VlnPlot(combined, group.by = 'region',features="nCount_RNA", pt.size=0, layer="counts") + NoLegend()
p2 = SpatialFeaturePlot(combined, features="nCount_RNA")
p1 | p2


#Look for bins with very few or very many counts. Bins with <30 counts or >2000 counts are suspicious and flagged for removal. 

### Bins flagged here should appear randomly across the tissue #####
#If they form a distinct pattern, reconsider the threshold — you might be removing biologically relevant data.

filtered_bins_counts = WhichCells(combined, expression=nCount_RNA < 25 | nCount_RNA > 2000)
SpatialDimPlot(combined, cells.highlight=filtered_bins_counts, cols.highlight=c("#FFFF00", "grey50")) + NoLegend() + ggtitle(paste(length(filtered_bins_counts), "bins filtered"))


p1 = VlnPlot(combined, features="nFeature_RNA", pt.size=0, layer="counts") + NoLegend()
p2 = SpatialFeaturePlot(combined, features="nFeature_RNA")
p1 | p2

filtered_bins_genes = WhichCells(combined, expression = nFeature_RNA < 20 | nFeature_RNA > 1000)
SpatialDimPlot(combined, cells.highlight=filtered_bins_genes, cols.highlight=c("#FFFF00", "grey50")) + NoLegend() + ggtitle(paste(length(filtered_bins_genes), "bins filtered"))


combined = PercentageFeatureSet(combined, pattern="^MT-", col.name="pMito_RNA") 

p1 = VlnPlot(combined, features="pMito_RNA", pt.size=0, layer="counts") + NoLegend()
p2 = SpatialFeaturePlot(combined, features="pMito_RNA") 
p1 | p2

filtered_bins_mito = WhichCells(combined, expression=pMito_RNA>20)
SpatialDimPlot(combined, cells.highlight=filtered_bins_mito, cols.highlight=c("#FFFF00", "grey50")) + NoLegend() + ggtitle(paste(length(filtered_bins_mito), "bins filtered"))


saveRDS(combined, file='Go_SpatialObject_simple.rds')

combined<-readRDS('Go_SpatialObject_simple.rds')

## Filtering of bins
#18085 104695  
dim(combined)
table(combined$region)

combined = subset(combined, nCount_RNA < 30 | nCount_RNA > 2000, invert=TRUE) %>% suppressWarnings()

#18085 87089
dim(combined)
table(combined$region)
#colon 49364 SI 37725 

saveRDS(combined, file='Go_SpatialObject_simple_clean.rds')

#combined<-readRDS('Go_SpatialObject_simple_clean.rds')

Idents(combined)<- combined$region
Colon<-subset(combined, idents='Colon')
SI<-subset(combined, idents='SI')

saveRDS(Colon, file='Go_SpatialObject_Colon_simple.rds')
saveRDS(SI, file='Go_SpatialObject_SI_simple.rds')


dim(Colon)
dim(SI)

## 3. Processing pipeline
combined<-SI
combined = NormalizeData(combined, normalization.method="LogNormalize", scale.factor=10000)
combined = FindVariableFeatures(combined, nfeatures=1000)
top10 = VariableFeatures(combined) %>% head(10)
top10

p = VariableFeaturePlot(combined, selection.method="vst")
p = LabelPoints(plot=p, points=top10, repel=TRUE)
p

combined = ScaleData(combined, features=VariableFeatures(combined))

combined = RunPCA(combined, reduction.name="pca", features=VariableFeatures(combined), npcs=80, nfeatures.print=5)
Reductions(combined, slot="pca")

DimPlot(combined, reduction="pca", dims=c(1, 2))
ElbowPlot(combined, reduction="pca", ndims=80)
VizDimLoadings(combined, reduction="pca", dims=1:4, nfeatures=10, balanced=TRUE)

combined = FindNeighbors(combined, reduction="pca", dims=1:40, k.param=20)
combined <- FindClusters(combined, resolution = 0.3, algorithm = 4, random.seed = 42)
table(combined$seurat_clusters, combined$orig.ident)
combined = RunUMAP(combined, reduction="pca", reduction.name="umap", dims=1:40, return.model=TRUE)

# cluster tree
for (res in c(0.1, 0.2, 0.3, 0.4,0.5,0.6, 0.7)) {
  combined <- FindClusters(combined, resolution = res, algorithm = 4, random.seed = 42,
                               verbose = FALSE)
}
pdf("ClusterTree_Spatial_Go_SI.pdf", width = 8, height = 8)
clustree(combined, prefix = 'RNA_snn_res.')
dev.off()


p1 = DimPlot(combined, reduction="umap", label=TRUE) + NoLegend()
p2 = SpatialDimPlot(combined) + NoLegend()
p1 | p2

#colon res=0.3
saveRDS(combined, file='Go_SpatialObject_SI_simple_clean_analyzed.rds')

#combined<-readRDS('Go_SpatialObject_simple_clean_analyzed.rds')


# Define the cluster to highlight
target_cluster <- "2:Paneth cells"

# Create a new variable that marks target vs others
combined_clean$highlight <- ifelse(combined_clean$Celltype == target_cluster,
                             paste0("Cluster_", target_cluster),
                             "Other")

# Plot
SpatialDimPlot(
  combined_clean,pt.size.factor = 2,
  group.by = "highlight",
  cols = c("Cluster_2:Paneth cells" = "red", "Other" = "gray80")
) +
  ggtitle(paste("Highlight Cluster", target_cluster))



##Cell calling
combined<-JoinLayers(combined)
markers = FindAllMarkers(object=combined,
                         test.use="wilcox",
                         layer="data",
                         only.pos=FALSE,
                         max.cells.per.ident=1000)

markers_top = markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=5) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

gt::gt(markers_top) %>% 
  gt::tab_options(container.height=450)



DotPlot(combined, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)")


VlnPlot(combined, features = c('nCount_RNA', 'nFeature_RNA', 'pMito_RNA'), pt.size = 0)

##General cell call colon
Idents(combined)<-combined$seurat_clusters
combined_clean<-subset(combined, idents=c('5', '13'), invert=TRUE)
ncells<-colnames(combined_clean)
combined<-readRDS('Go_SpatialObject_SI_simple.rds')
combined_clean<-subset(x = combined, cells = ncells)


combined_clean = NormalizeData(combined_clean, normalization.method="LogNormalize", scale.factor=10000)
combined_clean = FindVariableFeatures(combined_clean, nfeatures=1000)
combined_clean = ScaleData(combined_clean, features=VariableFeatures(combined_clean))

combined_clean = RunPCA(combined_clean, reduction.name="pca", features=VariableFeatures(combined_clean), npcs=80, nfeatures.print=5)
Reductions(combined_clean, slot="pca")

DimPlot(combined_clean, reduction="pca", dims=c(1, 2))
ElbowPlot(combined_clean, reduction="pca", ndims=80)
VizDimLoadings(combined_clean, reduction="pca", dims=1:4, nfeatures=10, balanced=TRUE)

combined_clean = FindNeighbors(combined_clean, reduction="pca", dims=1:60, k.param=20)
combined_clean <- FindClusters(combined_clean, resolution = 0.3, algorithm = 4, random.seed = 42)
combined_clean = RunUMAP(combined_clean, reduction="pca", reduction.name="umap", dims=1:40, return.model=TRUE)

# cluster tree
for (res in c(0.1, 0.2, 0.3, 0.4,0.5,0.6, 0.7)) {
  combined_clean <- FindClusters(combined_clean, resolution = res, algorithm = 4, random.seed = 42,
                           verbose = FALSE)
}
 

p1 = DimPlot(combined_clean, reduction="umap", label=TRUE) + NoLegend()
p2 = SpatialDimPlot(combined_clean) + NoLegend()
p1 | p2

p1
#colon res=0.4 and SI res=0.3
saveRDS(combined_clean, file='Go_SpatialObject_SI_simple_clean2_analyzed.rds')

combined_clean<-readRDS('Go_SpatialObject_SI_simple_clean2_analyzed.rds')
##Cell calling
combined_clean<-JoinLayers(combined_clean)
markers = FindAllMarkers(object=combined_clean,
                         test.use="wilcox",
                         layer="data",
                         only.pos=FALSE,
                         max.cells.per.ident=1000)

markers_top = markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=5) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

write.csv(markers_top, file='markers_top5_res03_SI_clean2.csv')

gt::gt(markers_top) %>% 
  gt::tab_options(container.height=450)



DotPlot(combined_clean, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)")


VlnPlot(combined_clean, features = c('nCount_RNA', 'nFeature_RNA', 'pMito_RNA'), pt.size = 0)

#Colon
Idents(combined_clean)<- combined_clean$seurat_clusters
new.cluster.ids <- c(
  '1:Fibroblast-like',
  '2:Enterocytes',
  '2:Enterocytes',
  '4:PCM',
  '5:Inflammatory macrophages',
  '6:Lymphocytes',
  '7:IgG Plasma cells',
  '6:Lymphocytes',
  '9:Enteroendocrine cells',
  '10:IgM/IgA Plasma cells',
  '11:Smooth muscle cells',
  '12:Goblet cells',
  '13:Crypt base progenitors',
  '14:BEST4 Enterocytes',
  '15:Inflammatory fibroblasts',
  '16:Smooth muscle / myofibroblast',
  '17:Inflammatory fibroblasts',
  '18:Endothelial (capillary)',
  '19:Endothelial (vascular)',
  '2:Enterocytes',
  '13:Crypt base progenitors',
  '22:Stromal fibroblasts',
  '23:Regenerative epithelium',
  '2:Enterocytes',
  '25:Activated fibroblast',
  '26:Proliferating epithelial progenitor',
  '7:IgG Plasma cells',
  '28:Stem cells',
  '29:Mast cells',
  '30:Smooth muscle subset',
  '12:Goblet cells',
  '32:Myofibroblasts',
  '33:Pericytes',
  '7:IgG Plasma cells',
  '35:Secretory progenitor ',
  '12:Goblet cells',
  '12:Goblet cells',
  '38:Smooth muscle',
  '38:Smooth muscle',
  '2:Enterocytes'
)
names(new.cluster.ids) <- levels(combined_clean)
combined_clean <- RenameIdents(combined_clean, new.cluster.ids)
combined_clean$Celltype <- Idents(combined_clean)


##SI 
new.cluster.ids <- c(
  '1:B cells',
  '2:Paneth cells',
  '3:IgA/IgG Plasma cells',
  '4:Enterocytes',
  '5:Stem cells',
  '6:IgG Plasma cells',
  '7:Enterocytes',
  '8:B cells',
  '4:Enterocytes',
  '2:Paneth cells',
  '4:Enterocytes',
  '4:Enterocytes',
  '13:Ileal enterocytes',
  '13:Ileal enterocytes',
  '15:Distal Enterocytes',
  '16:Enteroendocrine cells',
  '17:Macrophages',
  '16:Enteroendocrine cells',
  '19:IgA Plasma cells',
  '4:Enterocytes',
  '21:Smooth muscle',
  '22:Goblet cells',
  '4:Enterocytes',
  '24:Endothelial cells',
  '22:Goblet cells',
  '26:Stromal / fibroblast',
  '21:Smooth muscle',
  '28:Smooth muscle / myofibroblast',
  '29:TA epithelial ',
  '24:Endothelial cells',
  '31:Mast cells',
  '19:IgA Plasma cells',
  '33:M2 macrophages',
  '34:Fibroblasts / stromal',
  '35:Progenitor crypt epithelial',
  '4:Enterocytes',
  '37:Secretory progenitors',
  '29:TA epithelial',
  '19:IgA Plasma cells',
  '16:Enteroendocrine cells'
)
names(new.cluster.ids) <- levels(combined_clean)
combined_clean <- RenameIdents(combined_clean, new.cluster.ids)
combined_clean$Celltype <- Idents(combined_clean)


DimPlot(combined_clean, reduction="umap", label=TRUE, label.size = 3) + NoLegend()

saveRDS(combined_clean, file='Go_SpatialObject_SI_simple_clean2_analyzed.rds')

##SI





Colon<-readRDS('Go_SpatialObject_Colon_simple_clean2_analyzed.rds')
SI<-readRDS('Go_SpatialObject_SI_simple_clean2_analyzed.rds')

combined<-merge(Colon, SI)
##Plot
#Cluster interest
combined$Celltype_region<-paste(combined$Celltype, combined$region, sep='_')

genes_of_interest <- c("PGC", "TFF2", "AQP5", "MUC6", "CXCL17", "GP2", "BPIFB1", "TCN1")
genes_of_interest <- c( "AQP5", "MUC6",  "BPIFB1")
genes_of_interest <- c("WFDC2", "FAM3D")


genes_of_interest[!genes_of_interest %in% rownames(combined)]

combined <- AddModuleScore(
  object = combined,
  features = list(genes_of_interest),
  name = "Inflare"
)

Idents(combined) <- combined$Celltype_region  # or whatever column holds the cluster info
levels(combined)
combined$comparison_group <- ifelse(Idents(combined) == "4:PCM_Colon", "2:Paneth cells_SI", "Other")
table(combined$comparison_group)

VlnPlot(
  combined,
  features = "Inflare1",
  group.by = "comparison_group",
  pt.size = 0,layer  = 'scale.data',log = TRUE,
  cols = c("#1f77b4", "#aaaaaa")) +
  ggtitle("Module score for Inflare") +
  ylab("Module Score") +
  xlab("") +
  theme_bw()

combined$comparison_group <- case_when(
  Idents(combined) == "4:PCM_Colon" ~ "4:PCM_Colon",
  Idents(combined) == "2:Paneth cells_SI" ~ "2:Paneth cells_SI",
  TRUE ~ "Other"
)
table(combined$comparison_group)

pdf('ViolinPlot_PanethOlivier_PCM_Paneth_other.pdf', width = 5.5, height = 5)
VlnPlot(
  combined,
  features = "Inflare1",
  group.by = "comparison_group",
  pt.size = 0.1,layer  = 'scale.data',
  cols = c("#1f77b4", "#2ca02c", "#aaaaaa")) +
  ggtitle("Module score for PCM") +
  ylab("Module Score") +
  xlab("") +
  theme_bw()
dev.off()

dim(Colon)
table(Colon$Celltype)

library(patchwork)

# Assuming your Seurat object is called 'combined'
# and you already have cluster identities assigned
Idents(combined) <- "Celltype_region"  # or "seurat_clusters" if not yet renamed

# Subset to PCM_Colon cluster
pcm_subset <- subset(combined, idents = "4:PCM_Colon")

# get coordinates for slice1
coords4 <- GetTissueCoordinates(pcm_subset, image = "slice1")

# check range for proper zoom
range(coords4$x)
range(coords4$y)

# plot: full-resolution image with cluster 4 highlighted
pdf('SliceB21_Cluster4_PCM_Colon_REG3A.pdf', width = 12, height = 14)
SpatialPlot(
  pcm_subset,
  images = "slice1",
  features = 'REG3A',
  cols = c("4:PCM_Colon" = "red", "Other" = "gray80"),
  pt.size.factor = 2,         # adjust for HD density; smaller for denser data
  crop = FALSE,               # show full histology area
  alpha = 1                   # keep spots fully opaque
) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
dev.off()

# Check number of spots/cells
print(table(Idents(pcm_subset)))

# Visualize the spatial expression for antimicrobial markers
genes_of_interest <- c("DEFA5", "DEFA6", "REG3A", "PLA2G2A")

library(Seurat)
library(ggplot2)

genes_of_interest <- c("DEFA5", "DEFA6", "REG3A", "PLA2G2A")

# Generate plots for each gene and apply coord_cartesian to zoom in
plots <- lapply(genes_of_interest, function(gene) {
  p <- SpatialFeaturePlot(
    pcm_subset,
    features = gene,
    pt.size.factor = 4,
    crop = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    image.alpha = 0.8
  )
  # Apply zoom after plot creation
  p + 
    coord_cartesian(xlim = c(10000, 30000), ylim = c(15000, 35000)) +
    ggtitle(gene) +
    theme(plot.title = element_text(size = 14, face = "bold"))
})

# Combine all four into one grid
library(patchwork)
wrap_plots(plots)


pcm_subset = FindVariableFeatures(pcm_subset, nfeatures=1000)
pcm_subset = ScaleData(pcm_subset, features=VariableFeatures(pcm_subset))

pcm_subset = RunPCA(pcm_subset, reduction.name="pca", features=VariableFeatures(pcm_subset), npcs=80, nfeatures.print=5)
Reductions(pcm_subset, slot="pca")

DimPlot(pcm_subset, reduction="pca", dims=c(1, 2))
ElbowPlot(pcm_subset, reduction="pca", ndims=80)
VizDimLoadings(pcm_subset, reduction="pca", dims=1:4, nfeatures=10, balanced=TRUE)

pcm_subset = FindNeighbors(pcm_subset, reduction="pca", dims=1:80, k.param=20)
pcm_subset <- FindClusters(pcm_subset, resolution = 0.3, algorithm = 4, random.seed = 42)
pcm_subset = RunUMAP(pcm_subset, reduction="pca", reduction.name="umap", dims=1:40, return.model=TRUE)

# cluster tree
for (res in c(0.1, 0.2, 0.3, 0.4,0.5,0.6, 0.7)) {
  pcm_subset <- FindClusters(pcm_subset, resolution = res, algorithm = 4, random.seed = 42,
                                 verbose = FALSE)
}
pdf("ClusterTree_Spatial_Go_pcm_subset.pdf", width = 8, height = 8)
clustree(pcm_subset, prefix = 'RNA_snn_res.')
dev.off()


p1 = DimPlot(pcm_subset, reduction="umap", label=TRUE) + NoLegend()
p2 = SpatialDimPlot(pcm_subset) + NoLegend()
p1 | p2

p1

genes_of_interest <- c("DEFA5", "DEFA6", "REG3A", "PLA2G2A")
genes_of_interest <- c("DEFA5", "DEFA6", "REG3A", "PLA2G2A")

pdf("ViolinPlot_Spatial_Go_4PCM_subset_REG3A_DEFA5_DEFA6_PLA2G2A.pdf", width = 8, height = 8)
VlnPlot(
  pcm_subset,
  features = genes_of_interest,
  pt.size = 0,
  ncol = 2)
dev.off()

##Cell calling
pcm_subset<-JoinLayers(pcm_subset)
markers = FindAllMarkers(object=pcm_subset,
                         test.use="wilcox",
                         layer="data",
                         only.pos=FALSE,
                         max.cells.per.ident=1000)

markers_top = markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=5) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

write.csv(markers_top, file='markers_top5_res03_pcm_subset_clean2.csv')

gt::gt(markers_top) %>% 
  gt::tab_options(container.height=450)



DotPlot(pcm_subset, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)")


VlnPlot(pcm_subset, features = c('nCount_RNA', 'nFeature_RNA', 'pMito_RNA'), pt.size = 0)

Marker_genes<-c('IGHM',"CD74","ACTA2","STMN1",'JCHAIN',"LYZ",'REG3A',"FABP2","OLFM4",'UCHL1' , "MUC3A",  'MKI67',"DCLK1","CHGB",'LYVE1','MADCAM1','APOA4','MCT1','TPM1')

DotPlot(pcm_subset, features=Marker_genes %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)")

### Genes of PCM (signature)
##Cell calling
combined_clean<-JoinLayers(Colon)
Idents(Colon)<-Colon$Celltype
markers = FindAllMarkers(object=Colon,
                         test.use="wilcox",
                         layer="data",
                         only.pos=FALSE,
                         max.cells.per.ident=1000)

markers_sig<-subset(markers, markers$p_val_adj<0.05)
markers_sig_PCM<-subset(markers, markers$cluster %in% '4:PCM')

Idents(SI)<-SI$Celltype
markers = FindAllMarkers(object=SI,
                         test.use="wilcox",
                         layer="data",
                         only.pos=FALSE,
                         max.cells.per.ident=1000)

markers_sig<-subset(markers, markers$p_val_adj<0.05)
markers_sig_PC_SI<-subset(markers, markers$cluster %in% '2:Paneth cells')


head(markers_sig_PCM)
head(markers_sig_PC_SI)

# Genes in common between both dataframes
common_genes <- as.data.frame(intersect(markers_sig_PCM$gene, markers_sig_PC_SI$gene))
common_genes$Comparison<-'Common Paneth SI and PCM colon'
common_genes$GeneID<-common_genes$`intersect(markers_sig_PCM$gene, markers_sig_PC_SI$gene)`
common_genes$`intersect(markers_sig_PCM$gene, markers_sig_PC_SI$gene)`=NULL

# Genes exclusive to markers_sig_PCM
exclusive_PCM <- as.data.frame(setdiff(markers_sig_PCM$gene, markers_sig_PC_SI$gene))
exclusive_PCM$Comparison<-'PCM colon'
exclusive_PCM$GeneID<-exclusive_PCM$`setdiff(markers_sig_PCM$gene, markers_sig_PC_SI$gene)`
exclusive_PCM$`setdiff(markers_sig_PCM$gene, markers_sig_PC_SI$gene)`=NULL

# Genes exclusive to markers_sig_PanethSI
exclusive_SI <- as.data.frame(setdiff(markers_sig_PC_SI$gene, markers_sig_PCM$gene ))
exclusive_SI$Comparison<-'Paneth SI'
exclusive_SI$GeneID<-exclusive_SI$`setdiff(markers_sig_PC_SI$gene, markers_sig_PCM$gene)`
exclusive_SI$`setdiff(markers_sig_PC_SI$gene, markers_sig_PCM$gene)`=NULL

Table<-rbind(common_genes,exclusive_PCM, exclusive_SI )

write.csv(Table, 'SupplementaryTableX_SignatureGenes_PanetSI_PCMcolon.csv')

markers_top = markers_sig_PCM %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=20) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

markers_top$gene


write.csv(markers_top, file='markers_top5_res03_SI_clean2.csv')

gt::gt(markers_top) %>% 
  gt::tab_options(container.height=450)



DotPlot(combined_clean, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)")







Marker_genes<-c('IGHM',"CD74","ACTA2","STMN1",'JCHAIN',"LYZ",'REG3A',"FABP2","OLFM4",'UCHL1' , "MUC3A",  'MKI67',"DCLK1","CHGB",'LYVE1','MADCAM1','APOA4','MCT1','TPM1')


DotPlot(combined, features=Marker_genes %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="bottom") + 
  guides(size=guide_legend(order=1, title="Pct expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("In-house markers per cluster (expression scaled)")


set.seed(42)
cluster_bins = Cells(sp) %>% 
  sample(10000)


DoHeatmap(sp, features=Marker_genes, 
          cells=cluster_bins,
          group.colors=cal_pal50, 
          label=FALSE) + 
  scale_fill_viridis() +
  guides(colour="none")
```


# Cell type annotation by deconvolution

Spatial transcriptomics data often require deconvolution to **resolve the cellular composition** of individual spatial units (bins). This section guides you through the process of estimating cell type compositions in Visium HD data using the RCTD method from the `spacexr` package.

## What is deconvolution?

In spatial transcriptomics, such as 10x Visium HD, each bin may contain zero, one, or multiple cells, depending on the resolution, cell size, and local cell density. As a result, we cannot assume that each bin corresponds to a single cell. To address this, deconvolution methods are used to estimate the proportion of different cell types within each bin.

These methods also often incorporate spatial information (i.e., neighborhood relationships between bins) to improve accuracy.

One such method is **Robust Cell Type Decomposition (RCTD)**, implemented in the [`spacexr`](https://github.com/dmcable/spacexr) package. RCTD learns gene expression profiles of cell types from a **single-cell reference dataset** and uses these profiles to deconvolve spatial bins into estimated cell type compositions.

For this tutorial, we use a preprocessed reference derived from the [Tabula Muris senis](https://tabula-muris-senis.sf.czbiohub.org) single-cell dataset of large intestine.

We now prepare the `SpatialRNA` object required by `spacexr`:
  
  ```{r}
#| label: annot_1

# Counts
counts = GetAssayData(sp, assay="sketch", layer="counts")

# Coordinates
bins_sketch = Cells(sp[["sketch"]])
coords = GetTissueCoordinates(sp)[bins_sketch, 1:2]

# Total counts
total_counts = colSums(counts)

# Create spacexr SpatialRNA object
query = spacexr::SpatialRNA(coords, counts, total_counts)
```

## Load and prepare the reference dataset

We begin by loading the [Tabula Muris Senis large intestine dataset], which has been saved as a Seurat object:
  
  ```{r}
#| label: annot_2

load("datasets/TM_Intestine_facs_Figure.rda")
TM_Int<-UpdateSeuratObject(TM_Int)
TM_Int
```

We use the `Celltype` column for cell type annotations:
  
  ```{r}
#| label: annot_3

TM_Int[[]] %>% dplyr::select(Celltype) %>% head(10)
table(TM_Int$Celltype)
```

We convert this Seurat object into a `spacexr` Reference object:
  
  ```{r}
#| label: annot_4

# Counts
counts = GetAssayData(TM_Int, assay="RNA", layer="counts")

# Cell type labels
cell_types = as.factor(TM_Int$Celltype)
levels(cell_types) = gsub("/", "-", levels(cell_types))
cell_types = droplevels(cell_types)

# Total counts
total_counts = TM_Int$nReads

# Create Reference object
ref = spacexr::Reference(counts, cell_types, total_counts, n_max_cells=500)

gc()
```

## Run RCTD

We are now ready to run RCTD. However, since this step is time-consuming, we provide the results precomputed and saved in the file `datasets/RCTD.Rds`. For completeness, here is how you would normally run the algorithm:
  
  ```{r}
#| label: annot_5

# This will take five hours
#
RCTD = create.RCTD(query, ref, max_cores=12)
RCTD = run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD, "datasets/RCTD.Rds")
gc()
```

We simply load the results:
  
  ```{r}
#| label: annot_6

RCTD = readRDS("datasets/RCTD.Rds")
RCTD_results = RCTD@results$results_df
head(RCTD_results)
```

The `results_df` contains the main output of the deconvolution. Key columns include:  
  * `spot_class`: classification of each bin ("singlet", "doublet_certain", "doublet_uncertain", "reject")  
* `first_type`: predicted identity of the first cell type  
* `second_type`: only populated for doublets

We now add these predictions to the bin metadata of the sp object:
  
  ```{r}
#| label: annot_7

# Add to sp object
RCTD_results = RCTD_results[, c("spot_class", "first_type", "second_type")]
colnames(RCTD_results) = c("doublet", "first_cell_type", "second_cell_type")
sp = AddMetaData(sp, metadata=RCTD_results)

# Replace NA (not available) with "unknown"
sp$doublet = as.character(sp$doublet)
sp$doublet[is.na(sp$doublet)] = "unknown"
sp$doublet = factor(sp$doublet)

sp$first_cell_type = as.character(sp$first_cell_type)
sp$first_cell_type[is.na(sp$first_cell_type)] = "Unknown"
sp$first_cell_type = factor(sp$first_cell_type)

sp$second_cell_type = as.character(sp$second_cell_type)
sp$second_cell_type[is.na(sp$second_cell_type)] = "Unknown"
sp$second_cell_type = factor(sp$second_cell_type)
```



## Visualize the results

We can now visualize cell type predictions. For example, here we highlight two predicted secretory cells (typically labeled as Enteroendocrine, and Tuft cells) on the spatial plot:
  :
  ```{r}
#| label: annot_8

Idents(sp) = "first_cell_type"

secreted_cells = c("Paneth cells", "Enterocytes")
secreted_cells_bins = CellsByIdentities(sp, idents=secreted_cells)
secreted_cells_bins[["NA"]] = NULL

SpatialDimPlot(sp, cells.highlight=secreted_cells_bins, facet.highlight=TRUE, cols.highlight=c("#FFFF00", "grey50")) + NoLegend()
```
## Plot the reference-based cell annotations.
```{r}
#| label: plot_final
#| fig-height: 5
#| fig-width: 12

p1 = DimPlot(sp, reduction="umap.sketch", label=TRUE, group.by = 'first_cell_type') + NoLegend()+ scale_fill_manual(values=cal_pal50)
p2 = SpatialDimPlot(sp, group.by = 'first_cell_type') + scale_fill_manual(values=cal_pal50)
p1 | p2

```
## Annotate clusters and domains

Using the in house cell types, we can now try to annotate our clusters and domains. 

We plot the distribution of predicted cell types for the `seurat_clusters`:
  
  ```{r}
#| label: annot_10

# Call clusters and substitute cell types
Idents(sp)<-sp$seurat_clusters
new.cluster.ids <- c( "Immune cells", #1
                      "Stromal cells", #2
                      "Stem TA cells", #3
                      "Plasma cells", #4
                      "Paneth cells", #5
                      "Immature enterocytes", #6
                      "Enterocytes", #7
                      "Enterocytes", #8
                      "Enterocytes", #9
                      "Paneth cells", #10
                      "Smooth muscle", #11
                      "Enterocytes", #12
                      "Immature B cells", #13
                      "Tuft cells", #14
                      "Smooth muscl", #15
                      "Enteroendocrine", #16
                      "Pericyte", #17
                      "Immune cells", #18
                      "Enterocytes", #19
                      "Mast cells", #20
                      "Stromal cells" #21
)
names(new.cluster.ids) <- levels(sp)
sp <- RenameIdents(sp, new.cluster.ids)
sp$Celltype<-Idents(sp)

```

## Plot the cell annotations. always take into account that cell-types are subjects to granularity. Take immune cells for example, they have markers of B cells and macrophages, which 
```{r}
#| label: plot_final
#| fig-height: 5
#| fig-width: 12

p1 = DimPlot(sp, reduction="umap.sketch", label=TRUE) + NoLegend()+ scale_fill_manual(values=cal_pal50)
p2 = SpatialDimPlot(sp) + scale_fill_manual(values=cal_pal50)
p1 | p2

```
## Analysing the full dataset if you have no memory problems so far!

Once we have finished analyzing the sketched dataset, we can project the learned cluster labels, PCA, and UMAP from the sketch assay onto the entire dataset. This allows us to transfer the insights gained from the subsampled dataset to the full dataset, which may have additional complexity.

We use the `ProjectData` function to carry out this projection. Here’s the code to project the data from the sketch assay to the full dataset:
  
  ```{r}
#| label: plot_final2

sp = ProjectData(sp,
                 sketched.assay="sketch", assay="RNA", 
                 sketched.reduction="pca.sketch", full.reduction="pca", 
                 umap.model="umap.sketch",
                 refdata=list(seurat_clusters="seurat_clusters.sketch", Celltype.full='Celltype'),
                 dims=1:30)

# Fix the category levels for bin metadata column Seurat_clusters.008um and update the bin identities
sp$seurat_clusters = factor(sp$seurat_clusters, levels=levels(sp$seurat_clusters.sketch))
Idents(sp) = "seurat_clusters"

# Rename UMAP from full.umap.sketch to umap
umap = sp[["full.umap.sketch"]]
Key(umap) = "umap"
sp[["umap"]] = umap
sp[["full.umap.sketch"]] = NULL

# Set default assay back to RNA.008um
DefaultAssay(sp) = "RNA"
gc()
```

The Seurat object now contains a full PCA (`pca`), a full UMAP (`umap`) and a full clustering (`seurat_clusters). We color the clusters of the full dataset on the UMAP as well as in spatial context:
  
  ```{r}
#| label: plot_final3

p1 = DimPlot(sp, reduction="umap", label=TRUE, group.by = 'Celltype.full') +  scale_fill_manual(values=cal_pal50)
p2 = SpatialDimPlot(sp, group.by = 'Celltype.full') + scale_fill_manual(values=cal_pal50)
p1 | p2
```


```{r}
#| label: loupe_export

# Only for the 8 um bins
bins_to_export = Cells(sp[["RNA.008um"]])

# Bin metadata
metadata_to_export = sp[[]][bins_to_export, ] %>%
  dplyr::select(orig.ident, seurat_clusters.008um) %>%
  tibble::rownames_to_column("barcode")

write.csv(metadata_to_export, "datasets/exported_bin_metadata.csv", row.names=FALSE)

# UMAP
umap_to_export = Embeddings(sp[["umap.008um"]])[bins_to_export, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Barcode")
write.csv(umap_to_export, "datasets/exported_umap.csv", row.names=FALSE)
```

Information about importing results into Loupe can be found in the [Loupe tutorial](https://www.10xgenomics.com/support/software/loupe-browser/latest/tutorials/introduction/lb-sc-interface-and-navigation).






# Spatial clustering and domains

So far, we have explored the **fundamentals of analyzing spatial (and single-cell) transcriptomic data** using Seurat. The workflow introduced above serves is a good starting point.

However, a key limitation of this approach is that it treats each spatial bin independently, without considering its physical context within the tissue. Spatial analysis methods aim to address this by incorporating neighborhood information—enhancing biological interpretation and robustness. Several steps in our basic workflow can be extended to include this spatial context.


### If there is time: spatial information in clustering
## What is BANKSY?

[BANKSY](https://github.com/prabhakarlab/Banksy) is a computational method specifically designed for spatial transcriptomics analysis. By integrating information from neighboring bins, BANKSY can reduce noise in gene expression data, distinguish between different cell types based on spatial context, and identify spatial domains. 

BANKSY represents each cell or bin not only by its own transcriptomic profile but also by the average expression of its local neighborhood. This dual representation captures both intrinsic cellular identity and microenvironmental context. The degree to which neighborhood information is incorporated is controlled by the lambda parameter:  
  
  * Lower values (e.g., **lambda = 0.2**) emphasize immediate neighbors — useful for **cell typing**  
  * Higher values (e.g., **lambda = 0.8**) incorporate broader spatial context — useful for **domain detection**  
  
  Although BANKSY is a standalone tool with its own framework, it can be integrated into Seurat workflows using the `RunBanksy` function from the `SeuratWrappers` package. This allows users to apply BANKSY directly to Seurat objects.

In summary, BANKSY leverages the principle that a cell’s identity is shaped not only by its own gene expression but also by its surrounding environment. By adjusting a single hyperparameter, it provides flexibility to shift focus from cellular identity to higher-order tissue organization.

## Improve clustering with BANKSY

BANKSY can serve as an alternative to the default clustering approach we used earlier. In this step, we apply BANKSY with a low lambda value (0.2) to focus on local neighborhood information, and set the number of neighbors (k_geom) to 15.

**Note:** This step can be time- and memory-intensive.
This code sets up your R session to run computations in parallel using 10 cores and allows large objects (up to 100 GB) to be shared between the main session and the parallel workers.
```{r}
library(future)
options(future.globals.maxSize = 100 * 1024^3)  # must come before plan()
plan("multisession", workers = 10)

```


```{r}
#| label: spatial_1

sp = RunBanksy(sp,
               assay="RNA", assay_name="BANKSY",
               features="variable",
               lambda=0.2, k_geom=15,
               verbose=TRUE)

# Set it as default
DefaultAssay(sp) = "BANKSY"

gc()
```

At this point, we have created a new assay `BANKSY`, containing the BANKSY-transformed data, and have set it as the default.

**⚠️ Important**: BANKSY fills the `scale.data` layer. Therefore do not call `ScaleData` on the BANKSY assay as this negates the effects of lambda. 

```{r}
#| label: spatial_2

sp
```

Next, we apply PCA, clustering, and UMAP using the BANKSY assay. We use the first 15 PCs and set a clustering resolution of 0.5, though these values can and should be adjusted depending on your dataset.

```{r}
#| label: spatial_3

sp = RunPCA(sp, reduction.name="banksy_pca", features=rownames(sp[["BANKSY"]]), npcs=40, nfeatures.print=5)
sp = FindNeighbors(sp, reduction="banksy_pca", dims=1:30)
sp = FindClusters(sp, resolution=0.5, algorithm="leiden", random.seed=42)
sp = RunUMAP(sp, reduction="banksy_pca", reduction.name="banksy_umap", dims=1:40, return.model=TRUE)

sp$banksy_clusters = sp$seurat_clusters
DefaultAssay(sp) = "RNA"
```

Here we compare the results of default **Seurat clustering versus BANKSY clustering**, first using UMAP plots:
  
  ```{r}
#| label: spatial_4

p1 = DimPlot(sp, reduction="umap", label=TRUE, group.by="seurat_clusters") + NoLegend() + ggtitle("Seurat clustering")
p2 = DimPlot(sp, reduction="banksy_umap", label=TRUE, group.by="banksy_clusters") + NoLegend() + ggtitle("BANKSY clustering")
p1 | p2
```

And as spatial plots:
  
  ```{r}
#| label: spatial_5

p1 = SpatialDimPlot(sp, group.by="seurat_clusters") + NoLegend() + ggtitle("Seurat") + scale_fill_manual(values=cal_pal50)
p2 = SpatialDimPlot(sp, group.by="banksy_clusters") + NoLegend() + ggtitle("BANKSY") + scale_fill_manual(values=cal_pal50)
p1 | p2
```

## Recommendations
BANKSY does not automatically produce better results. Always evaluate both clustering outputs in the context of your biological question and data quality.

## Identifying spatial domains with BANKSY

In addition to improving clustering, BANKSY can be tuned to detect larger-scale tissue structures, often referred to as spatial domains. By increasing the `lambda` parameter and the number of neighbors (`k_geom`), BANKSY incorporates more spatial context and captures broader patterns of gene expression across the tissue.

In the following example, we set `lambda=0.8` and `k_geom=50` to include more neighboring bins in the analysis.

**Note:** We first remove the BANKSY assay from the previous run to conserve memory. This step is optional but will save us some memory resources.

```{r}
#| label: spatial_6

DefaultAssay(sp) = "RNA"
sp = RunBanksy(sp,
               assay="RNA", assay_name="BANKSY",
               features="variable",
               lambda=0.8, k_geom=50,
               verbose=TRUE)

# Set it as default
DefaultAssay(sp) = "BANKSY"

gc()
```

As before, we perform PCA, clustering, and UMAP on the BANKSY assay. We use the first 10 PCs and a clustering resolution of 0.5, but these values should again be adjusted depending on your dataset:
  
  ```{r}
#| label: spatial_7

sp = RunPCA(sp, reduction.name="banksy_pca", features=rownames(sp[["BANKSY"]]), npcs=40, nfeatures.print=5)
sp = FindNeighbors(sp, reduction="banksy_pca", dims=1:40)
sp = FindClusters(sp, resolution=0.5, algorithm="leiden")
sp$banksy_domains.008um = sp$seurat_clusters
```

If your clustering results look similar to those from the previous BANKSY run (with `lambda=0.2`), this could indicate that the dominant tissue structures are already captured at a local level. We encourage you to systematically explore different parameter settings for `lambda`, `k_geom`, and `resolution` to assess their impact on clustering outcomes and spatial domain resolution.

Since this BANKSY run is also complete, we again remove the assay and PCA reduction to free up memory:
  
  ```{r}
#| label: spatial_8

# Remove BANKSY assay and pca since not used anymore and we need the memory.
DefaultAssay(sp) = "RNA"
sp[["BANKSY"]] = NULL
sp[["banksy_pca<"]] = NULL
gc()
```

## Recommended reading

If you want to learn more about BANKSY, we find these pages helpful: 
  
  * https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md  
* https://github.com/prabhakarlab/Banksy?tab=readme-ov-file  


# Session info

To enhance reproducibility and facilitate the sharing of your analysis, it is good practice to include information about the R session and the packages used:
  
  ```{r}
#| label: session_info

sessioninfo::session_info()
```

# Useful resources

- Best practices:
  - [Single-cell best practices](https://www.sc-best-practices.org/preamble.html)
- [Orchestrating Single-Cell Analysis with Bioconductor
](https://bioconductor.org/books/release/OSCA/)
- [Orchestrating Spatial Transcriptomics Analysis with Bioconductor
](https://lmweber.org/OSTA/)
- Seurat:
  - [Essential commands](https://satijalab.org/seurat/articles/essential_commands)
- [Basic Single-cell Analysis Vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [Seurat Visium Analysis Vignette](https://satijalab.org/seurat/articles/spatial_vignette)
- [Seurat Visium HD Analysis Vignette](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette)
- Tools:
  - [Spaceranger](https://www.10xgenomics.com/support/software/space-ranger/latest)
- [Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser/latest)
- [Spotsweeper](https://github.com/MicTott/SpotSweeper)
- [banksy](https://github.com/prabhakarlab/Banksy)




