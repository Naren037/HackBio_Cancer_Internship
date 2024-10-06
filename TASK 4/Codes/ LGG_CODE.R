# README - https://github.com/Naren037/HackBio_Cancer_Internship/blob/87c6fe99a5e60f629a08c066d9278e19913d6e13/TASK%204/Codes/README

# START OF CODE 

# Required packages
library(gplots)
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(edgeR)
library(org.Hs.eg.db)
library(ggplot2)
library(DESeq2)
library(factoextra)
library(cluster) 
library(caret)
library(pROC)
library(DALEX)


############## UNSUPERVISED MACHINE LEARNING ##################


# Data Mining using TCGABiolinks
query_LGG <-GDCquery(project = 'TCGA-LGG',
                         data.category  = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         experimental.strategy = "RNA-Seq",
                         workflow.type = "STAR - Counts",)

GDCdownload(query_LGG)
DGE_LGG <- GDCprepare(query_LGG, summarizedExperiment = T)
LGG_Dataset <- assay(DGE_LGG, 'unstranded')
LGG_Dataset <- as.data.frame(LGG_Dataset)

# Preprocessing 1 - Variance stabilizing Transformation (VST)
pre1 <- DESeqDataSetFromMatrix(countData = LGG_Dataset, 
                              colData = colData(DGE_LGG), 
                              design = ~ 1)  # Unsupervised analysis, no outcome variable
pre1 <- DESeq(pre1)
pre1_data <- assay(vst(pre1))

# Pre-processing 2 - Selecting high variance genes
pre2 <- data.frame(t(pre1_data))
SD = apply(pre2, 2, sd)
Top_predictions = order(SD, decreasing = T)[1:10000]
pre2 = pre2[, Top_predictions]
dim(pre2)
pre2_data <- pre2

# Preprocessing 3 - Removing zero values  
pre3 <- preProcess(pre2_data, method = 'nzv', uniqueCut = 15)
pre3_data <- predict(pre3, pre2_data)
dim(pre3_data)

# Preprocessing 4 - Centering
pre4 <- preProcess(pre3_data, method = 'center')
pre4_data <- predict(pre4, pre3_data)
dim(pre4_data)

# Preprocessing 5 - Removing highly correlated features
correlation <- preProcess(pre4_data, method = 'corr', cutoff = 0.5) #ADJUST LATER 
Processed_data <- predict(correlation, pre4_data)
dim(Processed_data)

# Dimensionality Reduction - PCA
pca_res <- PCA((Processed_data), graph = FALSE)
# Visualize the first two principal components
fviz_pca_ind(pca_res,
             axes = c(1, 2),
             geom.ind = "point", 
             palette = "jco", 
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             legend.title = "Groups")

# K-Means Clustering
k <- 3
set.seed(123)
kmeans_res <- kmeans((Processed_data), centers = k)

# Add cluster labels to PCA plot
fviz_pca_ind(pca_res,
             geom.ind = "point", 
             col.ind = factor(kmeans_res$cluster),  
             palette = "jco",
             addEllipses = TRUE,
             #ellipse.type = "confidence",
             ellipse.level = 0.85,
             legend.title = "Cluster",
             title = "k-Means Clusterring")

# Double-check using Hierarchical Clustering
dist_matrix <- dist(Processed_data)
hc_res <- hclust(dist_matrix, method = "ward.D2")
# Plot dendrogram
plot(hc_res, labels = FALSE, hang = -1, main = "Hierarchical clustering", xlab = "Features")
rect.hclust(hc_res, k = 3, border = 2:4)

# Retrieving mutation data
Mutations_LGG <- data.frame(rownames(Processed_data), DGE_LGG$paper_IDH.codel.subtype, DGE_LGG$paper_ATRX.status, DGE_LGG$paper_TERT.expression.status)
rownames(Mutations_LGG) = Mutations_LGG$rownames.Processed_data.
Mutations_LGG$rownames.Processed_data. = NULL

# Subset clusters
clusters <- kmeans_res$cluster
kMeans_clusters <- data.frame(rownames(Processed_data), clusters, Mutations_LGG)
kMeans_clusters$rownames.Processed_data. = NULL
Cluster_1 <- subset(kMeans_clusters, clusters == 1)
Cluster_2 <- subset(kMeans_clusters, clusters == 2)
Cluster_3 <- subset(kMeans_clusters, clusters == 3)

# Export cluster details
write.csv(Cluster_1, "Cluster1.csv")
write.csv(Cluster_2, "Cluster2.csv")
write.csv(Cluster_3, "Cluster3.csv")

# Frequency of mutations in each cluster
# For cluster 1
IDH_cluster1 <- as.data.frame(table(Cluster_1$DGE_LGG.paper_IDH.codel.subtype))  
ATRX_cluster1 <- as.data.frame(table(Cluster_1$DGE_LGG.paper_ATRX.status)) 
TERT_cluster1 <- as.data.frame(table(Cluster_1$DGE_LGG.paper_TERT.expression.status))  
# For cluster 2
IDH_cluster2 <- as.data.frame(table(Cluster_2$DGE_LGG.paper_IDH.codel.subtype))  
ATRX_cluster2 <- as.data.frame(table(Cluster_2$DGE_LGG.paper_ATRX.status)) 
TERT_cluster2 <- as.data.frame(table(Cluster_2$DGE_LGG.paper_TERT.expression.status))
# For cluster 3
IDH_cluster3 <- as.data.frame(table(Cluster_3$DGE_LGG.paper_IDH.codel.subtype))  
ATRX_cluster3 <- as.data.frame(table(Cluster_3$DGE_LGG.paper_ATRX.status)) 
TERT_cluster3 <- as.data.frame(table(Cluster_3$DGE_LGG.paper_TERT.expression.status))

# Combined for mutations 
# IDH
IDH_Cluster <- data.frame(IDH_cluster1, IDH_cluster2$Freq, IDH_cluster3$Freq)
colnames(IDH_Cluster) = c("Mutation", "Cluster 1", "Cluster 2", "Cluster 3")
rownames(IDH_Cluster) = IDH_Cluster$Mutation
IDH_Cluster$Mutation = NULL
# ATRX
ATRX_Cluster <- data.frame(ATRX_cluster1, ATRX_cluster2$Freq, ATRX_cluster3$Freq)
colnames(ATRX_Cluster) = c("Mutation", "Cluster 1", "Cluster 2", "Cluster 3")
rownames(ATRX_Cluster) = ATRX_Cluster$Mutation
ATRX_Cluster$Mutation = NULL
# TERT
TERT_Cluster <- data.frame(TERT_cluster1, TERT_cluster2$Freq, TERT_cluster3$Freq)
colnames(TERT_Cluster) = c("Mutation", "Cluster 1", "Cluster 2", "Cluster 3")
rownames(TERT_Cluster) = TERT_Cluster$Mutation
TERT_Cluster$Mutation = NULL

# Export frequency results
write.csv(IDH_Cluster, "IDH_Cluster.csv")
write.csv(ATRX_Cluster, "ATRX_Cluster.csv")
write.csv(TERT_Cluster, "TERT_Cluster.csv")

# Bar plot
# IDH
barplot(IDH_Cluster, 
        legend.text = T,
        main = "IDH CLUSTERS",
        xlab = "Cluster group",
        ylab = "Number of features",
        ylim = c(0, 250))
# ATRX
barplot(ATRX_Cluster, 
        legend.text = T,
        main = "ATRX CLUSTERS",
        xlab = "Cluster group",
        ylab = "Number of features",
        ylim = c(0, 250))
# TERT
barplot(TERT_Cluster, 
        legend.text = T,
        main = "TERT CLUSTERS",
        xlab = "Cluster group",
        ylab = "Number of features",
        ylim = c(0, 250))


############## DGE ANALYSIS AND FUNCTIONAL ENRICHMENT ##################


DGE_Dataset <- t(LGG_Dataset)

# Pre-processing 1 - Adjusting NA values 
DGE_Dataset[is.na(DGE_Dataset)] <- rowMeans(DGE_Dataset, na.rm = TRUE)
dim(DGE_Dataset)

# Preprocessing 2 - Removing zero values  
DGE_Dataset_pre2 <- preProcess(DGE_Dataset, method = 'nzv', uniqueCut = 15)
DGE_Dataset <- predict(DGE_Dataset_pre2, DGE_Dataset)
dim(DGE_Dataset)

# Pre-processing 2 - Normalization using Trimmed mean of M values -> Adjusts differences in sample variation and sequencing depths
DGE_Dataset_TMM <- DGEList(counts = DGE_Dataset, group = NULL)
DGE_Dataset_TMM <- calcNormFactors(DGE_Dataset_TMM)
DGE_Dataset <- cpm(DGE_Dataset, log = F)
dim(DGE_Dataset)

# Pre-processing 3 - Upper Quantile filter of genes -> Removes lowly expressed genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = t(DGE_Dataset),
  method = "quantile",
  qnt.cut =  0.75)
dim(dataFilt)
DGE_Preprocessed <- t(dataFilt)
dim(DGE_Preprocessed)

# Retrieving IDH status
IDH_LGG <- data.frame(DGE_LGG$paper_IDH.status, rownames(Processed_data))
rownames(IDH_LGG) = IDH_LGG$rownames.Processed_data.
IDH_LGG$rownames.Processed_data. = NULL

# Adding IDH status
DGE_Preprocessed <- data.frame(IDH_LGG, DGE_Preprocessed)

# Subset by IDH status
LGG_WT <- subset(DGE_Preprocessed, DGE_LGG.paper_IDH.status == "WT")
dim(LGG_WT)
LGG_WT$DGE_LGG.paper_IDH.status = NULL
LGG_MUT <- subset(DGE_Preprocessed, DGE_LGG.paper_IDH.status == "Mutant")
dim(LGG_MUT)
LGG_MUT$DGE_LGG.paper_IDH.status = NULL

# Differential Gene Expression Analysis using TCGAnalyze with extremely stringent cut-offs for Biomarker discovery
LGG_WT <- t(LGG_WT)
LGG_MUT <- t(LGG_MUT)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = LGG_WT,
  mat2 = LGG_MUT,
  Cond1type = "WT",
  Cond2type = "Mutant",
  fdr.cut = 0.0005,
  logFC.cut = 4,
  method = "glmLRT")
plotdata <- dataDEGs[,2:6] # for further plots

# Retrieving gene names
ensembl_ids <- rownames(dataDEGs)
ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids)
gene_symbols <- mapIds(
  x = org.Hs.eg.db,           
  keys = ensembl_ids_clean,        
  column = "SYMBOL",          
  keytype = "ENSEMBL",       
  multiVals = "first")
dataDEGs <- data.frame(gene_symbols, dataDEGs)
dataDEGs$gene_name = NULL
dataDEGs$gene_type = NULL
print(gene_symbols)

# Heatmap of DEGs
# Getting expression data of DEGs
ensembl_ids <- rownames(dataFilt)
ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids)
rownames(dataFilt) = ensembl_ids_clean
heat.data <- dataFilt[rownames(dataDEGs),]
heat.data <- t(heat.data)
heat.data <- data.frame(DGE_Preprocessed$DGE_LGG.paper_IDH.status, heat.data)

# Heatmap.2
col_annotation <- ifelse(heat.data$DGE_Preprocessed.DGE_LGG.paper_IDH.status == "WT", "blue", "red")
heatmap.2(as.matrix(heat.data[,2:35]),
          Rowv = T,
          Colv = T,
          dendrogram = "both",
          scale = "row",
          trace = "none",
          col = bluered(200),
          RowSideColors = col_annotation,
          cexRow = 0.37,
          cexCol = 0.5,
          keysize = 1,
          key.title = 'Expression',
          margins = c(7,8),
          main = 'Differentially Expressed Genes',
          ylab = 'Sample Type',
          xlab = 'Genes')
legend("topright",  
       legend = c("WT", "Mutant"),  
       col = c("blue", "red"),  
       pch = 15,  
       cex = 1.2, 
       bty = "y")

# Bubble plot - DGE Analysis
ggplot(dataDEGs, aes(x = logFC, y = gene_symbols, size = -log10(PValue))) +
  geom_point(alpha=0.5)

# Exporting DGE analysis findings
write.csv(dataDEGs, "LGG_DEGs.csv")

# Functional Enrichment 
# Gene Ontology (GO) and Pathway enrichment by DEGs list

Genelist <- dataDEGs$gene_symbols
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes WT Vs Mutant",
  RegulonList = Genelist)


# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  GOBPTab = head(ansEA$ResBP),
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = Genelist,
  nBar = 10
)

# Print top 5 GO Biological processes
for (i in 1:5) { print(ansEA[["ResBP"]][[i]])}

# END OF CODE 
