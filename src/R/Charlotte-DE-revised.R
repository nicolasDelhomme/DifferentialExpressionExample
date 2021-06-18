#' ---
#' title: "Charlotte's DE (mod'ed)"
#' author: "Nicolas Delhomme, mod'ed from Charlotte Soneson's"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Temp. fix
options(connectionObserver = NULL)

#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(here)
  library(iSEE)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
})

#' * Helpers
source(here("UPSCb-common/src/R/volcanoPlot.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")

#' ```{r, load, echo=FALSE, eval=TRUE}
load(here("data/airway/dds.rda"))
load(here("data/airway/dge.rda"))
#' ```

#' # Differential Expression
#' ## DESeq2
#' Very easy, just one command
ds_se <- DESeq2::DESeq(ds_se)

#' ### Assessment
#' Here, we plot the dispersion estimation to ensure that the assumption that 
#' most genes are not differentially expressed holds
DESeq2::plotDispEsts(ds_se)

#' ### Results
#' Here we extract the results
res <- DESeq2::results(ds_se)

#' Wait! What was our design?
design(ds_se)

#' We have multiple possible contrast, which is the one R returns by default?
colData(ds_se)[,c("dex","cell")]
levels(ds_se$dex)
levels(ds_se$cell)

#' By default the last on in the that list (which is sorted alphabetically)
resultsNames(ds_se)

#' So we are looking at the effect of treated vs untreated. 
head(res)

#' OK, but what was my baseline?
levels(ds_se$cell)

#' Alright, we are looking at the effect of treated vs untreated in N052611
mcols(res, use.names = TRUE)

#' Let's take a global look at the results
summary(res)
hist(res$pvalue)

#' We add a couple of extra columns that will be useful for the interactive
#' visualization later
res$log10BaseMean <- log10(res$baseMean)
res$mlog10PValue <- -log10(res$pvalue)
rowData(ds_se)$log10Dispersion <- log10(rowData(ds_se)$dispersion)
rowData(ds_se)$DESeq2_dex_trt_vs_untrt <- res

#' ### Insights
#' #### Changing the FDR
res.05 <- results(ds_se, alpha = 0.05)
table(res.05$padj < 0.05)

#' #### Changing the log2FC - looking for a difference of expression of at least 1 fold change
resLFC1 <- results(ds_se, lfcThreshold = 1)
summary(resLFC1)
table(resLFC1$padj < 0.1)

#' #### Following Schurch et al., RNA, 2016 recommandations

#' ### Visualisation
#' #### Assessment
#' A volcano plot (the naming is probably obvious) is very useful to check
#' whether the assumption that most genes are not DE holds.
#' 
#' The plot has a representation of density as a color gradient 
#' from sparse to dense (gray -> blue -> red -> yellow). Note that the light blue, 
#' confusingly shows DE genes.
#' 
#' Here clearly most genes are located at 0,0 (the yellow dot)
volcanoPlot(res)

#' #### individual genes
plotCounts(ds_se, gene = "ENSG00000000003.14", intgroup = "dex", 
           normalized = TRUE, transform = FALSE)

#' ## EdgeR
#' Let us define our model
names(dge)
design <- model.matrix(~ cell + dex, data = dge$samples)

#' Perform the independent filtering (based on the model, filter the genes that have no 
#' power to be detected as differentially expressed. This aims at removing genes so lowly expressed 
#' that they represent "noise")
keep <- edgeR::filterByExpr(dge, design)
dge <- dge[keep, ]

#' Then we estimate the dispersion
dge <- edgeR::estimateDisp(dge, design)

#' And we can take a look at the dispersion estimation (dispersion == biological coefficient of varation)
edgeR::plotBCV(dge)

#' The, we can fit the model (using a quasi-likelihood f test)
fit <- edgeR::glmQLFit(dge, design)

#' and test
qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))

#' Finally looking at the results
tt.all <- edgeR::topTags(qlf, n = nrow(dge), sort.by = "none") # all genes
hist(tt.all$table$PValue)
tt <- edgeR::topTags(qlf, n = nrow(dge), p.value = 0.1) # genes with adj.p<0.1
tt10 <- edgeR::topTags(qlf) # just the top 10 by default
tt10

#' ## Comparison (DESeq2 as rows,edgeR as columns)
shared <- intersect(rownames(res), tt.all$table$gene.id)
table(DESeq2 = res$padj[match(shared, rownames(res))] < 0.1, 
      edgeR = tt.all$table$FDR[match(shared, tt.all$table$gene.id)] < 0.1)

#' including NAs to look at the independent filtering
table(DESeq2 = res$padj[match(shared, rownames(res))] < 0.1, 
      edgeR = tt.all$table$FDR[match(shared, tt.all$table$gene.id)] < 0.1,
      useNA="always")


#' let us do a visual comparison
grid.newpage()
d.sel <- match(shared, rownames(res))
e.sel <- match(shared, tt.all$table$gene.id)
grid.draw(venn.diagram(list(DESeq2 = rownames(res)[d.sel][res$padj[d.sel] < 0.1 & !is.na(res$padj[d.sel])], 
                       edgeR = tt.all$table$gene.id[e.sel][tt.all$table$FDR[e.sel] < 0.1 & !is.na(tt.all$table$FDR[e.sel])]),
          NULL,fill=pal[1:2]))

#' edgeR seems a little more lenient than DESeq2
sum(res$padj[match(shared, rownames(res))] < 0.1,na.rm = TRUE)

#' Looking at the differences in rank
plot(rank(res$pvalue[match(shared, rownames(res))]), 
     rank(tt.all$table$PValue[match(shared, tt.all$table$gene.id)]), 
     cex = 0.1, xlab = "DESeq2", ylab = "edgeR")

#' Same as for DESeq, looking for a more general hypothesis, lfc at least 1 or greater 
treatres <- edgeR::glmTreat(fit, coef = ncol(design), lfc = 1)
tt.treat <- edgeR::topTags(treatres, n = nrow(dge), sort.by = "none")

#' ## Another look at multiple testing
#' We have that many genes with a p.value below 0.05
sum(res$pvalue < 0.05, na.rm = TRUE)
#' and that many genes not NA
sum(!is.na(res$pvalue))

#' Assuming that there is no effect, we could expect that many genes
round(sum(!is.na(res$pvalue)) * 0.05)

#' expected ratio of false positives in the set with p < .05
round(sum(!is.na(res$pvalue))*0.05 / sum(res$pvalue < 0.05, na.rm = TRUE), 2)

#' How many genes below 10% FDR 
sum(res$padj < 0.1, na.rm = TRUE)

#' ## Visualisation
DESeq2::plotMA(res, ylim = c(-5, 5))

edgeR::plotSmear(qlf, de.tags = tt$table$gene.id)

#' ### variance stabilising transformation
vsd <- vst(ds_se,blind=FALSE)
mat <- assay(vsd)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("cell", "dex")])
pheatmap(mat, annotation_col = df)

#' ## Exporting the results
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = gsub("\\.[0-9]+$", "", row.names(res)),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

resOrdered <- res[order(res$padj), ]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[seq_len(100), ]
write.table(cbind(id = rownames(resOrderedDF), resOrderedDF), 
            file = here("analysis/DE/results.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE)

#' # Session Info
#'  ```{r, session info, echo=FALSE}
#' sessionInfo()
#' ```
