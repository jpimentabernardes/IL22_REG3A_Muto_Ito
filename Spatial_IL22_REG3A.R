---
  title: "Spatial RNAseq "
output: html_document
date: "2025-09-04"
---
  
library(Seurat)     
library(SeuratObject)
library(patchwork)  
library(tidyverse)  
library(grid)    
library(viridis) 
library(SeuratWrappers) 
library(Matrix)
library(readr)
library(dplyr)
library(tidyr)
library(arrow)  
library(jsonlite)
library(OpenImageR)

setwd('')

# Load dataset and information of tissue location
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


B22 = Load10X_Spatial("input/B22_07388/outs/binned_outputs/square_008um/", assay="RNA")
B22$orig.ident = "B22"
Idents(B22) = "B22"
B22

regionB22<-read.csv('input/B22_07388/outs/binned_outputs/square_008um/region.csv')

rownames(regionB22) <- regionB22$Barcode
regionB22 <- regionB22[colnames(B22), , drop = FALSE]
B22$region <- regionB22$region

# Merge the two slides
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

## check expression of genes of interest
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







