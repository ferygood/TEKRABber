#' Gene/TE expression data from human/chimpanzee brain RNA-seq
#' 
#' @description
#' Dataset contains 4 expression data from human and chimpanzee brain RNA-seq . 
#' We select raw fastq data from 10 humans and 10 chimpanzees from 
#' (Khrameeva E et al., 2020). Gene expression is generated using HISAT2 and 
#' featureCounts (Kim D et al., 2019; Liao Y et al., 2014). Transposable 
#' elements (TEs) expression is generated with multi-mapping option using STAR 
#' and TEtranscripts (Dobin A et al., 2013; Jin Y et al., 2015).
#' 
#' @format An object contains 4 expression counts:
#' \describe{
#'   \item{\strong{hmGene} human gene expression data}
#'   \item{\strong{hmTE} human TE expression}
#'   \item{\strong{chimpGene} chimpanzee gene expression data}
#'   \item{\strong{chimpTE} chimpanzee TE expression data}
#' }
#' 
#' @usage
#' data(speciesCounts)
#' hmGene <- speciesCounts$hmGene
#' hmTE <- speciesCounts$hmTE
#' chimpGene <- speciesCounts$chimpGene
#' chimpTE <- speciesCounts$chimpTE
"speciesCounts"


#' A subsets of normalized Gene/TE expression data from human/chimpanzee brain RNA-seq for 
#' correlation analysis demonstration
#' 
#' @description 
#' An object of class "TekCorrSet" which contains 4 expression counts. 
#' These data are generated from speciesCounts using TEKRABber. For a quick demo, 
#' we only select 50 orthologs and 50 transposable elements.
#' 
#' @format An object of class "TekCorrSet" which contains 4 expression counts 
#' and you can access it specifying the parameters using assay_tekcorrset():
#' \describe{
#'   \item{assay_tekcorrset(speciesCorr, "gene", "human")}{human gene expression data}
#'   \item{assay_tekcorrset(speciesCorr, "te", "human")}{human TE expression data}
#'   \item{assay_tekcorrset(speciesCorr, "gene", "chimpanzee")}{chimpanzee gene expression data}
#'   \item{assay_tekcorrset(speciesCorr, "te", "chimpanzee")}{chimpanzee TE expression data}
#' }
#' 
#' @usage 
#' data(speciesCorr)
#' hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
#' hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
#' chimpGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "chimpanzee")
#' chimpTECorrInput <- assay_tekcorrset(speciesCorr, "te", "chimpanzee")
"speciesCorr"

#' Input expression data of gene/TE for differentially expressed analysis within 
#' same species
#' 
#' @description 
#' TEKRABber can also be used comparing orthologs and transposable elements within
#' same species, i.e. control and treatment. Here we provide an example data for 
#' demonstration. This data based on syn8466812 RNA-seq (Allen M et al., 2016). 
#' However, the expression data is modified due to confidential agreement. 
#' Therefore, it can not represent the original data.
#' 
#' @format An object contains 2 expression data:
#' \describe{
#'   \item{\strong{gene}}{input gene data for DE analysis comparing control and treatment}
#'   \item{\strong{te}}{input TE data for DE analysis comparing control and treatment}
#' }
#' 
#' @usage 
#' data(ctInputDE)
#' geneInputDE <- ctInputDE$gene
#' teInputDE <- ctInputDE$te
#' 
"ctInputDE"

#' Normalized Gene/TE expression toy data in control and treatment in same species 
#' for correlation analysis
#' 
#' @description 
#' Dataset contains gene/TE expression data from control and treatment based on 
#' syn8466812 RNA-seq (Allen M et al., 2016) for correlation analysis. These data
#' are also modified due to confidential agreement. Therefore, it can not represent
#' the original data. For a quick demonstration, we only use 10 genes and 10 
#' transposable elements .
#' 
#' @format An object of class "TekCorrSet" which contains 4 expression counts 
#' and you can access it specifying the parameters using assay_tekcorrset():
#' \describe{
#'   \item{assay_tekcorrset(ctCorr, "gene", "control")}{control gene expression data}
#'   \item{assay_tekcorrset(ctCorr, "te", "control")}{control TE expression data}
#'   \item{assay_tekcorrset(ctCorr, "gene", "treatment")}{treatment gene expression data}
#'   \item{assay_tekcorrset(ctCorr, "te", "treatment")}{treatment TE expression data}
#' }
#' 
#' @usage 
#' data(ctCorr)
#' geneConCorrInput <- assay_tekcorrset(ctCorr, "gene", "control")
#' teConCorrInput <- assay_tekcorrset(ctCorr, "te", "control")
#' geneTreatCorrInput <- assay_tekcorrset(ctCorr, "gene", "treatment")
#' teTreatCorrInput <- assay_tekcorrset(ctCorr, "te", "treatment")
#' 
"ctCorr"

#' Orthology information of human and chimpanzee with scaling factor
#' 
#' @description 
#' An output list of data contains 2 elements after using orthologScale(). 
#' The first one is the orthology table comparing human and chimpanzee. The second one
#' is the scaling factor. The purpose of providing this dataset is to 
#' save time for user running the tutorial and give a template for demonstration.
#' 
#' @format An object contains 2 elements:
#' \describe{
#'   \item{\strong{orthologTable}}{containing orthology information retrieving from Ensembl}
#'   \item{\strong{scaleFactor}}{containing the scaling factor to normalize data}
#' }
#' 
#' @usage 
#' data(fetchDataHmChimp)
#' fetchData <- fetchDataHmChimp
#' fetchData$orthologTable
#' fetchData$scaleFactor
"fetchDataHmChimp"