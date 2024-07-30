.onLoad <- function(libname, pkgname) {
    messages <- c(
        "Welcome to TEKRABber version 1.8.0",
        "+ New function: prepareRMSK() for getting repeatmasker",
        "+ New parameter: `numCore` in corrOrthologTE() for parallel computing"
    )
    
    packageStartupMessage(paste(messages, collapse = "\n"))
}
