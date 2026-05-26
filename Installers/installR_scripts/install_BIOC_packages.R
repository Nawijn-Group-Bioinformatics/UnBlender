options(repos = c(
    CRAN = "https://cloud.r-project.org/"
))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install()

BiocManager::install(
    c("GenomeInfoDb", "SingleCellExperiment", "TOAST"),
    update=FALSE,
    ask=FALSE)
