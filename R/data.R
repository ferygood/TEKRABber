#' Gene/TE expression data from human/chimpanzee brain RNA-seq
#' 
#' Datasets containing 4 expression data from human and chimpanzee RNA-seq 
#' data. It is generated using in-house script using fastq data provided from
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
#' correlation analysis.
#' 
#' Datasets containing 2 normalized data from human and chimpanzee RNA-seq 
#' data. For a quick demo, here we only subset a 
#' 
#' @format A data list with 2 expression counts:
#' \describe{
#'   \item{geneCorr}{containing gene expression data from human and chimpanzee}
#'   \item{teCorr}{containing TE expression data from human and chimpanzee}
#' }
"speciesCorr"

#' Input expression data of gene/TE for differentially expressed analysis within same species
#' 
#' Datasets containing 2 expression toy data for demonstrate how to use TEKRABber
#' on experimetal design within the same species.
#' 
#' @format A data list with 2 expression counts:
#' \describe{
#'   \item{gene}{input gene data for DE analysis comparing control and treatment}
#'   \item{te}{input TE data for DE analysis comparing control and treatment}
#' }
"ctInputDE"

#' Normalized Gene/TE expression toy data in control and treatment in same species for 
#' correlation analysis.
#' 
#' Datasets containing 4 normalized toy data. For a quick demo, here we only subset a 
#' 
#' @format A data list with 4 expression counts:
#' \describe{
#'   \item{geneCorr}{containing gene expression data from control and treatment}
#'   \item{teCorr}{containing TE expression data from control and treatment}
#' }
"ctCorr"