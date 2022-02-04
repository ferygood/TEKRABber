library(TEKRABber)
library(shiny)
library(ggpubr)

DEresult = hmchimpDE
corrRef = hmCorrResult
corrCompare = chimpCorrResult
metadata = meta

geneInput <- DEresult$normalized_gene_counts
geneRes <- DEresult$gene_res
teInput <- DEresult$normalized_te_counts
teRes <- DEresult$te_res

rownames(geneInput) <- gsub("[()]", "", rownames(geneInput))
rownames(teInput) <- gsub("[()]", "", rownames(teInput))
rownames(geneRes) <- gsub("[()]", "", rownames(geneRes))
rownames(teRes) <- gsub("[()]", "", rownames(teRes))
corrRef$geneName <- gsub("[()]", "", corrRef$geneName)
corrRef$teName <- gsub("[()]", "", corrRef$teName)
corrCompare$geneName <- gsub("[()]", "", corrCompare$geneName)
corrCompare$teName <- gsub("[()]", "", corrCompare$teName)

# for expression tab: avoid subscript out of bounds with empty expression
geneChoices <- intersect(corrRef$geneName, rownames(geneInput))
geneChoices <- intersect(corrCompare$geneName, geneChoices)
teChoices <- intersect(corrRef$teName, rownames(teInput))
teChoices <- intersect(corrCompare$teName, teChoices)