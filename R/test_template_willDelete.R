

data("speciesCounts")
data("hg38_panTro6_rmsk")

hmGene <- speciesCounts$hmGene
hmTE <- speciesCounts$hmTE
chimpGene <- speciesCounts$chimpGene
chimpTE <- speciesCounts$chimpTE

fetchData <- orthologScale(
    speciesRef = "hsapiens",
    speciesCompare = "ptroglodytes",
    geneCountRef = hmGene,
    geneCountCompare = chimpGene,
    teCountRef = hmTE,
    teCountCompare = chimpTE,
    rmsk = hg38_panTro6_rmsk
)


inputBundle <- DECorrInputs(fetchData)

meta_data <- data.frame(
    species = c(rep("human", ncol(hmGene) -1), rep("chimpanzee", ncol(chimpGene) -1))
)

meta_data$species <- factor(meta_data$species, levels=c("human", "chimpanzee"))

rownames(meta_data) <- colnames(inputBundle$geneInputDESeq2)

testDE <- DEgeneTE(
    geneTable = inputBundle$geneInputDESeq2,
    teTable = inputBundle$teInputDESeq2,
    metadata = meta_data,
    expDesign = TRUE
)

#select 20 rows
hmGeneCorrInput <- testDE$geneCorrInputRef[c(1:20),]
hmTECorrInput <- testDE$teCorrInputRef[c(1:20),]

hmCorrResult <- corrOrthologTE(
    geneInput = hmGeneCorrInput,
    teInput = hmTECorrInput,
    corrMethod = "pearson",
    padjMethod = "fdr"
)

