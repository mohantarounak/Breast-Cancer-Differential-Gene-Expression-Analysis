if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "pheatmap", 
                       "EnhancedVolcano", "clusterProfiler", 
                       "org.Hs.eg.db", "AnnotationDbi"))

install.packages(c("ggplot2", "dplyr"))

library(GEOquery)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

gset <- getGEO(filename = "C:/Users/mohan/Desktop/bioinformatics project/GSE42568_series_matrix.txt/GSE42568_series_matrix.txt")
expr <- exprs(gset)
dim(expr)
head(expr)
dim(expr)

View(colon)
pdata <- pData(gset)
table(pdata$`characteristics_ch1`)

group <- c(rep("Normal", 17), rep("Tumor", 104))
group <- factor(group, levels = c("Normal", "Tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

write.csv(results, "DEG_results.csv", row.names = TRUE)
library(EnhancedVolcano)
EnhancedVolcano(results,
+                 lab = rownames(results),
+                 x = "logFC",
+                 y = "adj.P.Val",
+                 title = "Breast Cancer vs Normal",
+                 pCutoff = 0.05,
+                 FCcutoff = 1.5,
+                 pointSize = 2.0,
+                 labSize = 3.5)

library(pheatmap)
top50 <- rownames(results)[1:50]
expr_top50 <- expr[top50, ]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_top50)
 
pheatmap(expr_top50,
+          annotation_col = annotation_col,
+          show_rownames = FALSE,
+          show_colnames = FALSE,
+          scale = "row",
+          main = "Top 50 DEGs Heatmap",
+          color = colorRampPalette(c("blue","white","red"))(100))

BiocManager::install("hgu133plus2.db", ask = FALSE)
library(hgu133plus2.db)

probe_ids <- rownames(sig_genes)
gene_ids <- AnnotationDbi::select(hgu133plus2.db,
+                                 keys = probe_ids,
+                                 columns = c("SYMBOL", "ENTREZID"),
+                                 keytype = "PROBEID")

gene_ids <- gene_ids[!is.na(gene_ids$ENTREZID), ]
gene_ids <- gene_ids[!duplicated(gene_ids$ENTREZID), ]
head(gene_ids)

go_result <- enrichGO(gene = gene_ids$ENTREZID,
+                       OrgDb = org.Hs.eg.db,
+                       ont = "BP",
+                       pAdjustMethod = "BH",
+                       pvalueCutoff = 0.05)
dotplot(go_result, title = "GO Biological Process")
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,
+                           organism = "hsa",
+                           pvalueCutoff = 0.05)
go_result <- enrichGO(gene = gene_ids$ENTREZID,
+                       OrgDb = org.Hs.eg.db,
+                       ont = "BP",
+                       pAdjustMethod = "BH",
+                       pvalueCutoff = 0.05)
dotplot(go_result, title = "GO Biological Process")
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,
+                           organism = "hsa",
+                           pvalueCutoff = 0.05)
dotplot(kegg_result, title = "KEGG Pathway Enrichment")


png("GO_plot.png", width = 800, height = 600)
dotplot(go_result, title = "GO Biological Process")
dev.off()

png("KEGG_plot.png", width = 800, height = 600)
dotplot(kegg_result, title = "KEGG Pathway Enrichment")
dev.off()

png("Heatmap.png", width = 800, height = 600)
pheatmap(expr_top50,
+          annotation_col = annotation_col,
+          show_rownames = FALSE,
+          show_colnames = FALSE,
+          scale = "row",
+          main = "Top 50 DEGs Heatmap",
+          color = colorRampPalette(c("blue","white","red"))(100))
dev.off()

png("Volcano_plot.png", width = 800, height = 600)
EnhancedVolcano(results,
+                 lab = rownames(results),
+                 x = "logFC",
+                 y = "adj.P.Val",
+                 title = "Breast Cancer vs Normal",
+                 pCutoff = 0.05,
+                 FCcutoff = 1.5)
dev.off()

write.csv(results, "DEG_results.csv", row.names = TRUE)
write.csv(sig_genes, "significant_genes.csv", row.names = TRUE)
write.csv(gene_ids, "gene_ids.csv", row.names = TRUE)
list.files(getwd())