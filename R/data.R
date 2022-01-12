#' Gene/TE expression data from human/chimpanzee brain RNA-seq
#' 
#' Dataset contains 4 expression data from human and chimpanzee RNA-seq . 
#' It is generated using in-house script on fastq data provided from
#' (Khrameeva E et al., 2020).
#' 
#' @format A data list with 4 expression counts:
#' \describe{
#'   \item{hmGene}{human gene expression from RNA-seq}
#'   \item{hmTE}{human TE expression from RNA-seq}
#'   \item{chimpGene}{chimpanzee gene expression from RNA-seq}
#'   \item{chimpTE}{chimpanzee TE expression from RNA-seq}
#' }
"speciesCounts"

#' Normalized Gene/TE expression data from human/chimpanzee brain RNA-seq for 
#' correlation analysis
#' 
#' Dataset contains 2 normalized data from human and chimpanzee RNA-seq 
#' data. For a quick demo, we only use a subset of rows.
#' 
#' @format A data list with 2 expression counts:
#' \describe{
#'   \item{geneCorr}{containing gene expression data from human and chimpanzee}
#'   \item{teCorr}{containing TE expression data from human and chimpanzee}
#' }
"speciesCorr"

#' Input expression data of gene/TE for differentially expressed analysis within same species
#' 
#' Dataset contains 2 expression toy data for demonstrating how to use TEKRABber
#' on experimental design within the case of comparing same species.
#' 
#' @format A data list with 2 expression counts:
#' \describe{
#'   \item{gene}{input gene data for DE analysis comparing control and treatment}
#'   \item{te}{input TE data for DE analysis comparing control and treatment}
#' }
"ctInputDE"

#' Normalized Gene/TE expression toy data in control and treatment in same species for 
#' correlation analysis
#' 
#' Dataset contains 4 normalized toy data. For a quick demo, we only use a subset of rows.
#' 
#' @format A data list with 4 expression counts:
#' \describe{
#'   \item{geneCorr}{containing gene expression data from control and treatment}
#'   \item{teCorr}{containing TE expression data from control and treatment}
#' }
"ctCorr"

#' Orthology information of human and chimpanzee with scaling factor
#' 
#' A output list of data contains 2 elements after using orthologScale(). 
#' The first one is the orthology table comparing human and chimpanzee. The second one
#' is the scaling factor. The purpose of providing this dataset is to 
#' save time for user running the tutorial and give a template for demonstration.
#' 
#' @format A data list with 2 elements:
#' \describe{
#'   \item{orthologTable}{containing orthology information retrieving from Ensembl}
#'   \item{scaleFactor}{containing the scaling factor to normalize data}
#' }
"fetchDataHmChimp"