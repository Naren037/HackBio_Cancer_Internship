# READ.ME <- https://github.com/Naren037/hackbio-cancer-internship/blob/main/README_Task2.md

# START OF CODE

library(gplots) # Loading gplots to use heatmap.2
library(ggplot2) # Loading ggplot2 to make bubble plots

glioblastoma <- read.csv("C:/Users/Naren M/Desktop/Hackbio/Task 2/glioblastoma.csv", row.names=1) # Loading the data set

# "Generating heat maps"

# Diverging color palette - column clusterred 
heatmap.2(as.matrix(glioblastoma), Rowv = F, Colv = T, dendrogram = "col", scale = "row", trace = "none", col=hcl.colors(100, palette = 'Tofino'), margins = c(4,8), keysize = 2, key.title = 'Expression', lhei = c(0.8,5), lwid = c(1,2), cexRow = 0.67, cexCol = 0.96, xlab = 'SAMPLE', ylab = 'GENES', main = 'GLIOBLASTOMA')
# Sequential color palette - column clusterred 
heatmap.2(as.matrix(glioblastoma), Rowv = F, Colv = T, dendrogram = "col", scale = "row", trace = "none", col=hcl.colors(100, palette = 'Purples 3'), margins = c(4,8), keysize = 2, key.title = 'Expression', lhei = c(0.8,5), lwid = c(1,2), cexRow = 0.67, cexCol = 0.96, xlab = 'SAMPLE', ylab = 'GENES', main = 'GLIOBLASTOMA')
# Diverging color palette - row clusterred
heatmap.2(as.matrix(glioblastoma), Rowv = T, Colv = F, dendrogram = "row", scale = "row", trace = "none", col=hcl.colors(100, palette = 'Tofino'), margins = c(4,8), keysize = 2, key.title = 'Expression', lhei = c(0.8,5), lwid = c(1,2), cexRow = 0.67, cexCol = 0.96, xlab = 'SAMPLE', ylab = 'GENES', main = 'GLIOBLASTOMA')
# Diverging color palette - both clustered 
heatmap.2(as.matrix(glioblastoma), Rowv = T, Colv = T, dendrogram = "both", scale = "row", trace = "none", col=hcl.colors(100, palette = 'Tofino'), margins = c(4,8), keysize = 2, key.title = 'Expression', lhei = c(0.8,5), lwid = c(1,2), cexRow = 0.67, cexCol = 0.96, xlab = 'SAMPLE', ylab = 'GENES', main = 'GLIOBLASTOMA')

# "Grouping the samples based on heatmap cluster patterns"

group1 <- data.frame(glioblastoma$S1, glioblastoma$S2, glioblastoma$S3, glioblastoma$S4, glioblastoma$S5)
group2 <- data.frame(glioblastoma$S6, glioblastoma$S7, glioblastoma$S8, glioblastoma$S9, glioblastoma$S10) 

# "Foldchange calculation"

# Mean calculation
avg1 <- as.data.frame(rowMeans(group1))
avg2 <- as.data.frame(rowMeans(group2))
# Adjusting mean values of zero to elminate infinite fold change values
for (i in 1:nrow(avg1)) {if(avg1[i, ] == 0) {avg1[i, ] <- 1}}
for (j in 1:nrow(avg2)) {if(avg2[j, ] == 0) {avg2[j, ] <- 1}}  
Foldchange <- log2(avg1) - log2(avg2)

# "p-values calculation"

p_values <- numeric(nrow(glioblastoma))
for (i in 1:nrow(glioblastoma))
{      
  g1 <- as.numeric(group1[i, ])  
  g2 <- as.numeric(group2[i, ])  
  
  t_test_result <- t.test(g1, g2, var.equal = FALSE)
  
  p_values[i] <- t_test_result$p.value
}
p_values <- as.data.frame(p_values)

# "Finding differentially regulated genes"

genes <- data.frame(rownames(glioblastoma))
genes <- data.frame(genes, p_values, Foldchange)
upregulated_genes <- subset(genes, Foldchange > 1.5 & p_values < 0.05)
downregulated_genes <- subset(genes, Foldchange < -1.5 & p_values < 0.05)

# "Printing gene IDs"

downregulated_genes$rownames.glioblastoma
upregulated_genes$rownames.glioblastoma

# "Bubble plot"

ggplot(enrichment_all, aes(x = log10.p.value, y = Pathway, size = No.of.Genes)) +
  +     geom_point(alpha=0.5)

# END OF CODE
