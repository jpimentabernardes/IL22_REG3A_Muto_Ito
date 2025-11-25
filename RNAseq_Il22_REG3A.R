library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(variancePartition)
library(RColorBrewer)
library(pheatmap)  
library(clusterProfiler)  
library(org.Hs.eg.db)  
library(forcats)

# set window
setwd('')

# set seed
set.seed(42)

# Load quantification files (Salmon output)
files <- list.files(path = "output_files", pattern = "_quant.sf", full.names = TRUE)

# Load sample metadata
metadata <- read.csv('/Users/joana/Documents/IKMB/Manuscripts/Go_Revision_NatureCommu/Metadata_IL22_REG3A.csv', stringsAsFactors = FALSE)

# Load tx2gene mapping and clean it (remove second column)
tx2gene <- read.table("output_files/tx2gene.tsv")
tx2gene$V2 <- NULL

# Extract sample names from filenames, e.g. "/path/REG3A1_quant.sf" -> "REG3A1"
sample_names <- gsub(".*/|_quant.sf", "", files)

# Match files order to metadata sample order
files <- files[match(metadata$sample, sample_names)]

# Import transcript counts using tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)



# Create DESeq dataset with no design (~1)
dds <- DESeqDataSetFromTximport(txi, metadata, ~1)

# Filter genes: keep genes with non-zero counts in at least 20% of samples
dds <- dds[rowSums(counts(dds) == 0) < 0.2 * ncol(txi$counts), ]

# Estimate size factors (normalization)
dds <- estimateSizeFactors(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Rename columns with sample names
colnames(normalized_counts) <- metadata$sample

# Variance Stabilizing Transformation
vst_counts <- vst(dds)
vst_mat <- assay(vst_counts)

# Rename columns for consistency
colnames(vst_mat) <- metadata$sample

# Principal Component Analysis (PCA)
pc <- prcomp(t(vst_mat))
pc_df <- as.data.frame(pc$x)
pc_df <- cbind(pc_df, metadata)

txi_counts<-as.data.frame(txi$counts)
write.csv(txi_counts, file = 'output2/Raw_counts.csv')
write.csv(vst_mat, file = 'output2/Normalized_counts.csv')

# Summary of PCA for axis labels
pc_summary <- summary(pc)
table(pc_df$organoid)

# Plot PCA (PC1 vs PC3) colored by Stimulation
pdf("output2/02_pca_plot_colored_by_PC12_Stimulation.pdf", width = 7, height = 5)
P <- ggplot(data = pc_df, mapping = aes(x = PC1, y = PC2, fill = Stimulation, shape = organoid)) +
  geom_point(size = 5,color = "black") +
  xlab(paste0("PC1 (", round(pc_summary$importance["Proportion of Variance", "PC1"] * 100, 2), "% variance)")) +
  ylab(paste0("PC2 (", round(pc_summary$importance["Proportion of Variance", "PC2"] * 100, 2), "% variance)")) +
  theme_bw() +
  scale_fill_manual(
    name = "Code",
    values = c("CTRL" = "#7083a4", "IL22" = "#f95b70")) +
  scale_shape_manual(
    name = "organoid",
    values = c("SF11" = 21, "SF13" = 24, 'SY8'=22)) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22))
print(P)
dev.off()



#Plotting residuals instead
library(limma)

normalized_counts <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ 1) %>%
  estimateSizeFactors() %>%
  counts(normalized = TRUE) %>%
  as.data.frame()

metadata$Stimulation
#remove batch effect of sequencing run for visualization
fit <- lmFit(normalized_counts, model.matrix(~ Stimulation, metadata))
normalized_counts_residuals <- residuals(fit, normalized_counts)
colnames(normalized_counts_residuals) <- metadata$sample

# Genes of interest for heatmap
genes_of_interest <- c("CH507-254M2.2", "FAM72C", "REG3A", "ARFGAP1", 
                       "CASZ1", "LAMB2", "LINC00649", "MUC12", 
                       "MUC5B", "SART1", "SIPA1", 'MUC6')

# Subset vst matrix for genes of interest
subset_vst <- normalized_counts_residuals[rownames(normalized_counts_residuals) %in% genes_of_interest, ]

# Prepare annotation for heatmap: only the Stimulation column
annotation_col <- metadata
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[, "Stimulation", drop = FALSE]

#rownames(metadata) <- rownames(annotation_col) <- colnames(subset_vst)

# Order your samples by Stimulation
order_cols <- order(metadata$Stimulation)

# Reorder both your expression matrix and the annotation data frame
subset_vst <- subset_vst[, order_cols]
annotation_col <- as.data.frame(annotation_col[order_cols, , drop = FALSE])

# Double-check rownames match
stopifnot(all(rownames(annotation_col) == colnames(subset_vst)))

# Plot heatmap of selected genes with annotation
pdf('output2/DEG_Stimulation_Go_Rev_REG3A_new_MUC6_residuals.pdf')
pheatmap(subset_vst, 
  scale = "row", 
  annotation_col = annotation_col,
  cluster_rows = TRUE,           # keep row order
  cluster_cols = FALSE,           # keep column order
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Heatmap of Selected Genes"
)
dev.off()


###Delta change
library(dplyr)
library(stringr)

# Step 1 – Extract replicate number from the sample name
metadata <- metadata %>%
  mutate(replicate = str_extract(sample, "_\\d+$") %>% str_remove("_"))

# Step 2 – Verify column order matches metadata
stopifnot(all(colnames(normalized_counts) %in% metadata$sample))

# Reorder metadata to match normalized_counts
metadata <- metadata[match(colnames(normalized_counts), metadata$sample), ]

# Step 3 – Split data by organoid and replicate
# For each unique (organoid, replicate) pair, find the CTRL and IL22 samples
delta_list <- list()

for (org in unique(metadata$organoid)) {
  for (rep in unique(metadata$replicate)) {
    
    # Subset metadata for this pair
    sub_meta <- metadata %>%
      filter(organoid == org, replicate == rep)
    
    # Skip if we don't have both CTRL and IL22 samples
    if (!all(c("CTRL", "IL22") %in% sub_meta$Stimulation)) next
    
    # Get sample names
    ctrl_sample <- sub_meta$sample[sub_meta$Stimulation == "CTRL"]
    il22_sample <- sub_meta$sample[sub_meta$Stimulation == "IL22"]
    
    # Compute delta for each gene
    delta <- normalized_counts[, il22_sample] - normalized_counts[, ctrl_sample]
    
    # Name the column clearly, e.g., "SY8_rep1"
    delta_list[[paste0(org, "_rep", rep)]] <- delta
  }
}

# Step 4 – Combine all results into one data frame
delta_df <- as.data.frame(delta_list)
rownames(delta_df) <- rownames(normalized_counts)

head(delta_df)

# Genes of interest for heatmap
genes_of_interest <- c("CH507-254M2.2", "FAM72C", "REG3A", "ARFGAP1", 
                       "CASZ1", "LAMB2", "LINC00649", "MUC12", 
                       "MUC5B", "SART1", "SIPA1", 'MUC6')

genes_of_interest <- c('DEFA5','DEFA6','MUC6')
genes_of_interest <- c( "AQP5", "MUC6",  "BPIFB1")


# Subset vst matrix for genes of interest
subset_vst <- delta_df[rownames(delta_df) %in% genes_of_interest, ]

pdf('output2/DEG_Stimulation_Go_Rev_REG3A_new_MUC6_DeltaChange.pdf')
pheatmap(subset_vst, 
         scale = "row", 
        #annotation_col = annotation_col,
         cluster_rows = TRUE,           # keep row order
         cluster_cols = FALSE,           # keep column order
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Heatmap of Selected Genes"
)
dev.off()





significant<-c('PRSS2', 'DMBT1', 'REG1A', 'CFI', 'PLA2G2A',
      'MATN2', 'IL1R2', 'INSL5', 'GFRA3', 'TOX2')


# Subset vst matrix for genes of interest
subset_vst <- vst_mat[rownames(vst_mat) %in% significant, ]

# Prepare annotation for heatmap: only the Stimulation column
annotation_col <- metadata
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[, "Stimulation", drop = FALSE]

# Reorder annotation rows to match subset_vst columns exactly
annotation_col <- annotation_col[colnames(subset_vst), , drop = FALSE]

# Check if names match (this should pass)
stopifnot(all(colnames(subset_vst) == rownames(annotation_col)))

# Plot heatmap of selected genes with annotation
pdf('output2/DEG_Top5_Go_Rev_REG3A_new.pdf')
pheatmap(subset_vst, 
         scale = "row", 
         annotation_col = annotation_col,
         cluster_rows = FALSE,              # Disable row clustering
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Heatmap of Selected Genes")
dev.off()


Inflare<-c('REG3A','DEFA6','DEFA5', 'PLA2G2A', 'WFDC2', 'RPL18A', 'FAM83D', 'VSIG2', 'ICAM3', 'FXYD3', 'INTS1', 'ADAMTS6')

res_df$Gene<-rownames(res_df)
Inflare<-c('MUC6', 'PGC', 'AQP5',  'BPIFB1')
# Subset vst matrix for genes of interest
subset_vst <- vst_mat[rownames(vst_mat) %in% Inflare, ]

# Prepare annotation for heatmap: only the Stimulation column
annotation_col <- metadata
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[, "Stimulation", drop = FALSE]

# Reorder annotation rows to match subset_vst columns exactly
annotation_col <- annotation_col[colnames(subset_vst), , drop = FALSE]

# Check if names match (this should pass)
stopifnot(all(colnames(subset_vst) == rownames(annotation_col)))

# Plot heatmap of selected genes with annotation
pdf('output2/Inflare2_Stimulation_Go_Rev_REG3A_new.pdf')
pheatmap(subset_vst, 
         scale = "row", 
         annotation_col = annotation_col,
         cluster_rows = FALSE,        
         clustering_distance_rows = "euclidean",
         # Disable row clustering
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Heatmap of Selected Genes")
dev.off()

## Differential expression
dds <- DESeqDataSetFromTximport(txi,  metadata,  design = ~ organoid + Stimulation )  

# Apply variance stabilizing transformation if needed
dds <- DESeq(dds)

# Get results for Flare (1 vs. 0)
res <- results(dds)
res <- res[order(res$padj), ]
res_df <- as.data.frame(res)

significant<-subset(res_df, res_df$padj < 0.05)
significant$ID<-rownames(significant)

significant_Up<-subset(significant, significant$log2FoldChange > 0)
significant_Down<-subset(significant, significant$log2FoldChange < 0)
significant_Up$signal<-'Upregulated'
significant_Down$signal<-'Downregulated'
significant_Genes<-rbind(significant_Up, significant_Down)

write.csv(significant_Genes, file='output2/DEGs_Stimulated_controlvsIL22_organoids.csv')
present_genes_entrez <- bitr(rownames(res_df), 
                             fromType = "SYMBOL", 
                             toType="ENTREZID", 
                             OrgDb=org.Hs.eg.db)$ENTREZID

# gene ontology enrichment analysis based on overrepresentation
GOresults.list <- list()
for(i in unique(significant_Genes$signal)){
  print(i)
  markers <- significant_Genes[significant_Genes$signal==i,]
  print(paste("cluster: ",i, ", marker genes: ",paste(markers$ENTREZID[1:10],collapse=", "),", ...",sep=""))
  
  markers_entrez <- bitr(markers$ID, 
                         fromType = "SYMBOL", 
                         toType="ENTREZID", 
                         OrgDb=org.Hs.eg.db)$ENTREZID
  
  # GO enrichment
  GO <- as.data.frame(enrichGO(gene = markers_entrez,
                               universe = present_genes_entrez,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "bonferroni",
                               keyType = "ENTREZID",
                               pvalueCutoff  = 0.2,
                               qvalueCutoff  = 0.2,
                               readable      = T))
  
  if(nrow(GO)>0){GO$cluster <- as.character(i)}
  if(nrow(GO)>0){GOresults.list[[paste(i)]] <- GO}
}



#### GO
GO <- do.call("rbind", GOresults.list)

CompGOresults_sub<-subset(GO, GO$p.adjust <0.05)
CompGOresults_sub2<-subset(CompGOresults_sub, CompGOresults_sub$Count>=5)

write.csv(CompGOresults_sub2, file='output2/GO_DEG_Go_REG_rev_Stimulation_organoid.csv')


# Plot Go enerichment analysis
CompGOresults_sub2 <- CompGOresults_sub2 %>%
  mutate(Description = fct_reorder2(Description, cluster, qvalue, .desc = FALSE))

pdf("output2/GO_DEG_Go_REG_rev_Stimulation_organoid.pdf", width = 12, height = 4.5)
ggplot(CompGOresults_sub2, aes(x = qvalue, y = Description, fill=cluster)) +
  geom_point(aes(size = Count), shape=21) +
  #facet_grid( ~ cluster.x, , space = "free_y") +
  ylab(NULL) +
 scale_fill_manual(values = c(  'cornflowerblue',"brown1"))+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
