\name{TEKRABber}
\alias{TEKRABber}
\docType{package}
\title{
  An R package estimates the correlations of orthologs and transposable elements between two species
}
\description{
  TEKRABber is made to provide an user-friendly pipeline for comparing 
    orthologs and transposable elements (TEs) between two species. It considers the 
    orthology confidence between two species from BioMart to normalize expression 
    counts and detect differentially expressed orthologs/TEs. Then it provides one 
    to one correlation analysis for desired orthologs and TEs. There is also an app 
    function to have a first insight on the result. Users can prepare orthologs/TEs 
    RNA-seq expression data by their own preference to run TEKRABber following the 
    data structure mentioned in the vignettes.
}
\details{
  TEKRABber analysis pipeline includes 4 main functions:
  
  1. \strong{orthologScale()}: obtain orthology information and calculate scaling factor.
  
  2. \strong{DECorrInputs()}: create the input files for running DE / correlation analysis.
  
  3. \strong{DEgeneTE()}: run DE analysis on orthologs and transposable elements.
  
  4. \strong{corrOrthologTE()}: estimate correlation between selected orthologs and transposable elements.
  
  5. \strong{appTEKRABber()}: (optional) find first insight from data using an local webapp.
  
  Find more details in vignette or on the helping page, i.e. ?orthologScale
}

\author{
Yao-Chung Chen, Katja Nowick.

Maintainer: Yao-Chung Chen \email{yao-chung.chen@@fu-berlin.de}

\href{https://github.com/ferygood/TEKRABber}{TEKRABber GitHub Repo}

}