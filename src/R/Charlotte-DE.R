ds_se <- DESeq2::DESeq(ds_se)
DESeq2::plotDispEsts(ds_se)

res <- DESeq2::results(ds_se)
head(res)

mcols(res, use.names = TRUE)

summary(res)
hist(res$pvalue)

## We also add a couple of extra columns that will be useful for the interactive
## visualization later
res$log10BaseMean <- log10(res$baseMean)
res$mlog10PValue <- -log10(res$pvalue)
rowData(ds_se)$log10Dispersion <- log10(rowData(ds_se)$dispersion)
rowData(ds_se)$DESeq2_dex_trt_vs_untrt <- res

res.05 <- results(ds_se, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(ds_se, lfcThreshold = 1)
summary(resLFC1)
table(resLFC1$padj < 0.1)

plotCounts(ds_se, gene = "ENSG00000000003.14", intgroup = "dex", 
           normalized = TRUE, transform = FALSE)

names(dge)

design <- model.matrix(~ cell + dex, data = dge$samples)

keep <- edgeR::filterByExpr(dge, design)
dge <- dge[keep, ]
dge <- edgeR::estimateDisp(dge, design)
edgeR::plotBCV(dge)

fit <- edgeR::glmQLFit(dge, design)
qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
tt.all <- edgeR::topTags(qlf, n = nrow(dge), sort.by = "none") # all genes
hist(tt.all$table$PValue)
tt <- edgeR::topTags(qlf, n = nrow(dge), p.value = 0.1) # genes with adj.p<0.1
tt10 <- edgeR::topTags(qlf) # just the top 10 by default
tt10

shared <- intersect(rownames(res), tt.all$table$gene.id)
table(DESeq2 = res$padj[match(shared, rownames(res))] < 0.1, 
      edgeR = tt.all$table$FDR[match(shared, tt.all$table$gene.id)] < 0.1)

sum(res$padj[match(shared, rownames(res))] < 0.1,na.rm = TRUE)

plot(rank(res$pvalue[match(shared, rownames(res))]), 
     rank(tt.all$table$PValue[match(shared, tt.all$table$gene.id)]), 
     cex = 0.1, xlab = "DESeq2", ylab = "edgeR")

treatres <- edgeR::glmTreat(fit, coef = ncol(design), lfc = 1)
tt.treat <- edgeR::topTags(treatres, n = nrow(dge), sort.by = "none")

sum(res$pvalue < 0.05, na.rm = TRUE)
sum(!is.na(res$pvalue))

round(sum(!is.na(res$pvalue)) * 0.05)

round(sum(!is.na(res$pvalue)) * 0.05) # expected 'null' less than 0.05
sum(res$pvalue < 0.05, na.rm = TRUE) # observed p < .05
# expected ratio of false positives in the set with p < .05
round(sum(!is.na(res$pvalue))*0.05 / sum(res$pvalue < 0.05, na.rm = TRUE), 2)

sum(res$padj < 0.1, na.rm = TRUE)

DESeq2::plotMA(res, ylim = c(-5, 5))

edgeR::plotSmear(qlf, de.tags = tt$table$gene.id)

suppressPackageStartupMessages({
  library(pheatmap)
})
mat <- assay(vsd)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("cell", "dex")])
pheatmap(mat, annotation_col = df)

suppressPackageStartupMessages({
  library(iSEE)
})

#app <- iSEE(ds_se)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys = gsub("\\.[0-9]+$", "", row.names(res)),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

resOrdered <- res[order(res$padj), ]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[seq_len(100), ]
write.table(cbind(id = rownames(resOrderedDF), resOrderedDF), 
            file = "results.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)
