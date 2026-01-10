.onLoad <- function(libname, pkgname) {
    messages <- c(
        "Welcome to TEKRABber version 1.15.0",
        "+ prepareRMSK(): automatically retrieves RepeatMasker data via AnnotationHub",
        "+ corrOrthologTE(): new parameter `numCore` for parallel computing"
    )

    packageStartupMessage(paste(messages, collapse = "\n"))
}
