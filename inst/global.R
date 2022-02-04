library(TEKRABber)
library(shiny)
library(ggpubr)

geneInput <- appDE$normalized_gene_counts
geneRes <- appDE$gene_res
teInput <- appDE$normalized_te_counts
teRes <- appDE$te_res

rownames(geneInput) <- gsub("[()]", "", rownames(geneInput))
rownames(teInput) <- gsub("[()]", "", rownames(teInput))
rownames(geneRes) <- gsub("[()]", "", rownames(geneRes))
rownames(teRes) <- gsub("[()]", "", rownames(teRes))
appRef$geneName <- gsub("[()]", "", appRef$geneName)
appRef$teName <- gsub("[()]", "", appRef$teName)
appCompare$geneName <- gsub("[()]", "", corrRef$geneName)
appCompare$teName <- gsub("[()]", "", corrRef$teName)

# for expression tab: avoid subscript out of bounds with empty expression
geneChoices <- intersect(appRef$geneName, rownames(geneInput))
geneChoices <- intersect(corrRef$geneName, geneChoices)
teChoices <- intersect(appRef$teName, rownames(teInput))
teChoices <- intersect(corrRef$teName, teChoices)