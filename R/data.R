#' Gene/TE expression data from human/chimpanzee brain RNA-seq
#' 
#' @description
#' Dataset contains 4 expression data from human and chimpanzee brain RNA-seq. 
#' We select raw fastq data from 10 humans and 10 chimpanzees from 
#' (Khrameeva E et al., 2020). Gene expression is generated using HISAT2 and 
#' featureCounts (Kim D et al., 2019; Liao Y et al., 2014). Transposable 
#' elements (TEs) expression is generated with multi-mapping option using 
#' STAR and TEtranscripts (Dobin A et al., 2013; Jin Y et al., 2015).
#' @usage data(speciesCounts)
#' @format An object contains 4 expression counts:
#' \describe{
#'     \item{hmGene}{human gene expression data}
#'     \item{hmTE}{human TE expression}
#'     \item{chimpGene}{chimpanzee gene expression data}
#'     \item{chimpTE}{chimpanzee TE expression data}
#' }
#' 
#' @examples
#' data(speciesCounts)
#' hmGene <- speciesCounts$hmGene
#' hmTE <- speciesCounts$hmTE
#' chimpGene <- speciesCounts$chimpGene
#' chimpTE <- speciesCounts$chimpTE
"speciesCounts"

#' Input expression data of gene/TE for differentially expressed analysis 
#' within same species
#' 
#' @description 
#' TEKRABber can also be used comparing orthologs and transposable elements 
#' within same species, i.e., control and treatment. Here we provide an 
#' example data for demonstration. This data was based on syn8466812 
#' RNA-seq (Allen M et al., 2016). However, the expression data was modified 
#' due to confidential agreement. Therefore, it cannot represent 
#' the original data.
#' @usage data(ctInputDE)
#' @format An object contains 2 expression data:
#' \describe{
#'     \item{\strong{gene}}{
#'     input gene data for DE analysis comparing control and treatment}
#'     \item{\strong{te}}{
#'     input TE data for DE analysis comparing control and treatment}
#' }
#' 
#' @examples 
#' data(ctInputDE)
#' geneInputDE <- ctInputDE$gene
#' teInputDE <- ctInputDE$te
#' 
"ctInputDE"

#' Example output comparing human and chimpanzee data using orhtologScale()
#' 
#' @description 
#' An output list of data contains 7 elements after using orthologScale(), 
#' including (1) orthology table comparing human and chimpanzee. (2) scaling 
#' factor for orthologous genes (3) gene count table from reference species (4)
#'  gnee count table from species you want to compare (5) scaling factor for TEs
#' (6) TE count table from reference species (7) TE count table from the 
#' species you want to compare. The aim to provide this dataset is to save 
#' time for user running the vignettes and give a 
#' template for demonstration.
#' @usage data(fetchDataHmChimp)
#' @format An object contains 2 elements:
#' \describe{
#'     \item{\strong{orthologTable}}{orthology information from Ensembl}
#'     \item{\strong{scaleFactor}}{scaling factor to normalize data}
#' }
#' 
#' @examples 
#' data(fetchDataHmChimp)
#' fetchData <- fetchDataHmChimp
#' fetchData$orthologTable
#' fetchData$scaleFactor
"fetchDataHmChimp"

#' Repeatmasker track annotations with human and chimpanzee
#' 
#' @description 
#' This Repeatmasker track annotations table was first downloaded from UCSC 
#' Genome Table Browser and it included the name, class, and average gene 
#' length in repeats(transposable elements). This data is used for demonstrate 
#' an example for user how to provide a annotation table to normalize their 
#' data which in this case comparing human(hg38) to chimpanzee(panTro6).
#' @usage data(hg38_panTro6_rmsk)
#' @examples 
#' data(hg38_panTro6_rmsk)
"hg38_panTro6_rmsk"
