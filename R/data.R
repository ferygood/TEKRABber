#' Gene expression dataset from human brain RNA-seq data
#' 
#' A dataset containing randomly slected 10 human brain cells bulked RNA-seq
#' gene expression data. It is generated using in-house script by the author 
#' using fastq data provided from (Khrameeva E et al., 2020).
#' 
#' @format A data frame with 64252 rows and 11 variables:
#' \describe{
#'   \item{Geneid}{Ensembl gene ID}
#'   \item{SRR8750453}{sample ID}
#'   \item{SRR8750454}{sample ID}
#'   \item{SRR8750455}{sample ID}
#'   \item{SRR8750456}{sample ID}
#'   \item{SRR8750457}{sample ID}
#'   \item{SRR8750458}{sample ID}
#'   \item{SRR8750459}{sample ID}
#'   \item{SRR8750460}{sample ID}
#'   \item{SRR8750461}{sample ID}
#'   \item{SRR8750462}{sample ID}
#' }
"hmGene"

#' Transposable elements expression dataset from human brain RNA-seq data
#' 
#' A dataset containing randomly slected 10 human brain cells bulked RNA-seq
#' transposable elements expression data. It is generated using in-house script 
#' by the author using fastq data provided from (Khrameeva E et al., 2020).
#' 
#' @format A data frame with 992 rows and 11 variables:
#' \describe{
#'   \item{teName}{Name of transposable elements}
#'   \item{SRR8750453}{sample ID}
#'   \item{SRR8750454}{sample ID}
#'   \item{SRR8750455}{sample ID}
#'   \item{SRR8750456}{sample ID}
#'   \item{SRR8750457}{sample ID}
#'   \item{SRR8750458}{sample ID}
#'   \item{SRR8750459}{sample ID}
#'   \item{SRR8750460}{sample ID}
#'   \item{SRR8750461}{sample ID}
#'   \item{SRR8750462}{sample ID}
#' }
"hmTE"

#' Gene expression dataset from human brain RNA-seq data
#' 
#' A dataset containing randomly slected 10 chimpanzee brain cells bulked 
#' RNA-seq gene expression data. It is generated using in-house script by the 
#' author using fastq data provided from (Khrameeva E et al., 2020).
#' 
#' @format A data frame with 16922 rows and 11 variables:
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{SRR8750637}{sample ID}
#'   \item{SRR8750638}{sample ID}
#'   \item{SRR8750639}{sample ID}
#'   \item{SRR8750640}{sample ID}
#'   \item{SRR8750641}{sample ID}
#'   \item{SRR8750642}{sample ID}
#'   \item{SRR8750643}{sample ID}
#'   \item{SRR8750644}{sample ID}
#'   \item{SRR8750645}{sample ID}
#'   \item{SRR8750646}{sample ID}
#' }
"chimpGene"

#' Transposable elements expression dataset from chimpanzee brain RNA-seq data
#' 
#' A dataset containing randomly slected 10 chimpanzee brain cells bulked 
#' RNA-seq transposable elements expression data. It is generated using in-house
#' script by the author using fastq data provided from (Khrameeva E et al., 2020).
#' 
#' @format A data frame with 1260 rows and 11 variables:
#' \describe{
#'   \item{teName}{Name of transposable elements}
#'   \item{SRR8750637}{sample ID}
#'   \item{SRR8750638}{sample ID}
#'   \item{SRR8750639}{sample ID}
#'   \item{SRR8750640}{sample ID}
#'   \item{SRR8750641}{sample ID}
#'   \item{SRR8750642}{sample ID}
#'   \item{SRR8750643}{sample ID}
#'   \item{SRR8750644}{sample ID}
#'   \item{SRR8750645}{sample ID}
#'   \item{SRR8750646}{sample ID}
#' }
"chimpTE"

#' Human gene dataset for running correlation analysis
#' 
#' A normalized gene expression dataset containing 10 human samples to perform 
#' a quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 50 rows and 10 variables:
#' \describe{
#'   \item{SRR8750453}{sample ID}
#'   \item{SRR8750454}{sample ID}
#'   \item{SRR8750455}{sample ID}
#'   \item{SRR8750456}{sample ID}
#'   \item{SRR8750457}{sample ID}
#'   \item{SRR8750458}{sample ID}
#'   \item{SRR8750459}{sample ID}
#'   \item{SRR8750460}{sample ID}
#'   \item{SRR8750461}{sample ID}
#'   \item{SRR8750462}{sample ID}  
#' }
#' 
"hmGeneCorrInput"

#' Human transposable elements dataset for running correlation analysis
#' 
#' A normalized transposable elements expression dataset containing 10 human 
#' samples to perform a quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 50 rows and 10 variables:
#' \describe{
#'   \item{SRR8750453}{sample ID}
#'   \item{SRR8750454}{sample ID}
#'   \item{SRR8750455}{sample ID}
#'   \item{SRR8750456}{sample ID}
#'   \item{SRR8750457}{sample ID}
#'   \item{SRR8750458}{sample ID}
#'   \item{SRR8750459}{sample ID}
#'   \item{SRR8750460}{sample ID}
#'   \item{SRR8750461}{sample ID}
#'   \item{SRR8750462}{sample ID}  
#' }
#' 
"hmTECorrInput"

#' Chimpanzee gene dataset for running correlation analysis
#' 
#' A normalized gene expression dataset containing 10 chimpanzee samples to perform 
#' a quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 50 rows and 10 variables:
#' \describe{
#'   \item{SRR8750637}{sample ID}
#'   \item{SRR8750638}{sample ID}
#'   \item{SRR8750639}{sample ID}
#'   \item{SRR8750640}{sample ID}
#'   \item{SRR8750641}{sample ID}
#'   \item{SRR8750642}{sample ID}
#'   \item{SRR8750643}{sample ID}
#'   \item{SRR8750644}{sample ID}
#'   \item{SRR8750645}{sample ID}
#'   \item{SRR8750646}{sample ID}
#' }
#' 
"chimpGeneCorrInput"

#' Chimpanzee transposable elements dataset for running correlation analysis
#' 
#' A normalized transposable elements expression dataset containing 10 chimpanzee 
#' samples to perform a quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 50 rows and 10 variables:
#' \describe{
#'   \item{SRR8750637}{sample ID}
#'   \item{SRR8750638}{sample ID}
#'   \item{SRR8750639}{sample ID}
#'   \item{SRR8750640}{sample ID}
#'   \item{SRR8750641}{sample ID}
#'   \item{SRR8750642}{sample ID}
#'   \item{SRR8750643}{sample ID}
#'   \item{SRR8750644}{sample ID}
#'   \item{SRR8750645}{sample ID}
#'   \item{SRR8750646}{sample ID}  
#' }
#' 
"chimpTECorrInput"

#' Gene expression dataset for running differentially expressed analysis (DE)
#' 
#' A toy gene expression dataset containing 5 controls and 5 treatments to demonstrate
#' the case if user want to use TEKRABber for DE analysis in the same species.
#'
#' @format A data frame with 33223 rows and 10 variables:
#' \describe{
#'   \item{control_1}{control sample ID}
#'   \item{control_2}{control sample ID}
#'   \item{control_3}{control sample ID}
#'   \item{control_4}{control sample ID}
#'   \item{control_5}{control sample ID}
#'   \item{treatment_6}{treatment sample ID}
#'   \item{treatment_7}{treatment sample ID}
#'   \item{treatment_8}{treatment sample ID}
#'   \item{treatment_9}{treatment sample ID}
#'   \item{treatment_10}{treatment sample ID}
#' }
"geneInputDE"

#' Transposable elements expression dataset for running differentially expressed 
#' analysis (DE)
#' 
#' A toy transposable elements expression dataset containing 5 controls and 
#' 5 treatments to demonstrate the case if user want to use TEKRABber for DE 
#' analysis in the same species.
#'
#' @format A data frame with 687 rows and 10 variables:
#' \describe{
#'   \item{control_1}{control sample ID}
#'   \item{control_2}{control sample ID}
#'   \item{control_3}{control sample ID}
#'   \item{control_4}{control sample ID}
#'   \item{control_5}{control sample ID}
#'   \item{treatment_6}{treatment sample ID}
#'   \item{treatment_7}{treatment sample ID}
#'   \item{treatment_8}{treatment sample ID}
#'   \item{treatment_9}{treatment sample ID}
#'   \item{treatment_10}{treatment sample ID}
#' }
"teInputDE"

#' Gene expression data in control group for running correlation analysis
#' 
#' A small toy gene expression dataset containing 5 controls to perform quick 
#' demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 10 rows and 5 variables:
#' \describe{
#'   \item{control_1}{control sample ID}
#'   \item{control_2}{control sample ID}
#'   \item{control_3}{control sample ID}
#'   \item{control_4}{control sample ID}
#'   \item{control_5}{control sample ID}
#' }
"geneConCorrInput"

#' Transposable elements expression data in control group for running correlation analysis
#'
#' A small toy transposable elements expression dataset containing 5 controls to
#' perform quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 10 rows and 5 variables:
#' \describe{
#'   \item{control_1}{control sample ID}
#'   \item{control_2}{control sample ID}
#'   \item{control_3}{control sample ID}
#'   \item{control_4}{control sample ID}
#'   \item{control_5}{control sample ID}
#' }
"teConCorrInput"

#' Gene expression data in treatment group for running correlation analysis
#' 
#' A small toy gene expression dataset containing 5 treatments to perform quick 
#' demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 10 rows and 5 variables:
#' \describe{
#'   \item{treatment_6}{control sample ID}
#'   \item{treatment_7}{control sample ID}
#'   \item{treatment_8}{control sample ID}
#'   \item{treatment_9}{control sample ID}
#'   \item{treatment_10}{control sample ID}
#' }
"geneTreatCorrInput"

#' Transposable elements expression data in treatment group for running 
#' correlation analysis
#' 
#' A small toy transposable elements expression dataset containing 5 treatments 
#' to perform quick demonstration using TEKRABber for correlation analysis.
#' 
#' @format A data frame with 10 rows and 5 variables:
#' \describe{
#'   \item{treatment_6}{control sample ID}
#'   \item{treatment_7}{control sample ID}
#'   \item{treatment_8}{control sample ID}
#'   \item{treatment_9}{control sample ID}
#'   \item{treatment_10}{control sample ID}
#' }
"teTreatCorrInput"