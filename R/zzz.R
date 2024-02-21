.onLoad <- function(libname, pkgname) {
    messages <- c(
        "Welcome to TEKRABber version 1.8.0",
        "+ New function: prepareRMSK() for getting rmsk track",
        "+ Convert Ensembl ID to gene name in appTEKRABber()",
        "+ with setting number of cores for parallel computing"
    )
    
    packageStartupMessage(paste(messages, collapse = "\n"))
}
