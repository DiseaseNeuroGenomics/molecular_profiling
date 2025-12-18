---
title: "remaCor analysis"
output: html_document
---

```{r setup, message=FALSE}
#libraries and outdirs

library(edgeR)
library(variancePartition)
library(remaCor)
library(Matrix)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(data.table)
library(cowplot)
library(ggdendro)
library(limma)
library(reshape2)

dir.create("remacor/plots_1", recursive = TRUE, showWarnings = FALSE)
dir.create("remacor/plots_2", recursive = TRUE, showWarnings = FALSE)

```
```{r}
#function for running remacor
#requires voom object, fit object (from dream/limma),gene name,transcript info)
#runs per gene
get_mvTest = function(gene_name,
fit_obj,
voom_obj,
Coef,
TranscriptInfo) {

Gene_features = which(TranscriptInfo$ensembl_gene_id == gene_name)

CollinearTranscriptFlag = "No"
combTransName = TransKeep = TransRem = NA

trans_cor = cor(t(voom_obj$E[Gene_features,]))
diag(trans_cor) = 0
CollinearTranscripts =
rownames(trans_cor)[which(apply(trans_cor, 1, max) > 0.999)]

if (length(CollinearTranscripts) > 1) {
combTransName = paste(CollinearTranscripts, collapse = "_")
TransRem = CollinearTranscripts[-1]
Gene_features =
setdiff(Gene_features,
which(rownames(TranscriptInfo) %in% TransRem))
CollinearTranscriptFlag = "Yes"
}

mvTest_results =
c(
mvTest(fit = fit_obj, vobj = voom_obj,
features = Gene_features, coef = Coef,
method = "RE2C")[c("stat.FE", "stat.het", "pvalue")],
mvTest(fit = fit_obj, vobj = voom_obj,
features = Gene_features, coef = Coef,
method = "FE")[c("stat", "pvalue", "n_features")],
CollinearTranscriptFlag,
combTransName,
CollinearTranscripts[1],
paste(TransRem, collapse = "_")
)

names(mvTest_results) =
c("stat_mean_RE2C", "stat_hetero_RE2C", "pvalue_RE2C",
"stat_FE", "pvalue_FE", "transcript_count",
"CollinearTranscriptFlag",
"Collinear_Transcript_Names",
"Transcript_Kept",
"Transcript_Skipped")

mvTest_results
}
```


```{r}
### load and process data
load("~/placement/transcript_analysis/DET_Analysis_ALL.Rdata")

tx_df = data.frame(
modeledVoomObj$genes$PeakID,
modeledVoomObj$genes$transcript_name,
modeledVoomObj$genes$gene_name
)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

tx_ens_gene =
getBM(
attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
filters = "ensembl_transcript_id",
values = tx_df$modeledVoomObj.genes.PeakID,
mart = mart
)

id = match(tx_df$modeledVoomObj.genes.PeakID,
tx_ens_gene$ensembl_transcript_id)

tx_df2 = tx_df
tx_df2$ensembl_gene_id = tx_ens_gene$ensembl_gene_id[id]
rownames(tx_df2) = tx_df2$modeledVoomObj.genes.PeakID


```
```{r}
##run analysis
 	genes_to_test = unique(tx_df2$ensembl_gene_id)

results =
lapply(
na.omit(genes_to_test),
get_mvTest,
fit_obj = fitDream,
voom_obj = modeledVoomObj,
Coef = "celltype.SCZ_Control",
TranscriptInfo = tx_df2
)

get_resdf = function(x, rnames) {
x = rbindlist(x)
x$rnames = na.omit(rnames)
x = unique(x)
x = as.data.frame(x)
rownames(x) = x$rnames
x$rnames = NULL
x
}

res2 = get_resdf(celltype_results, genes_to_test)
saveRDS(res2, "remacor/celltype_results.RDS")


```



```{r}
##plot functions
plotCor_RKmod = function(cor) {
ID = rownames(cor)
df2 = melt(cor)
colnames(df2)[1:2] = c("Var1", "Var2")
df2$Var1 = factor(df2$Var1, ID)
df2$Var2 = factor(df2$Var2, rev(ID))
df2$label = round(df2$value, 1)

ggplot(df2, aes(Var1, Var2, fill = value)) +
geom_tile() +
geom_text(aes(label = label), size = 3) +
scale_fill_gradient2(
low = "blue", mid = "white", high = "red",
limits = c(-1, 1)
) +
theme_minimal() +
theme(
panel.grid = element_blank(),
axis.text.x = element_text(angle = 90, hjust = 1),
aspect.ratio = 1
)
}
```

