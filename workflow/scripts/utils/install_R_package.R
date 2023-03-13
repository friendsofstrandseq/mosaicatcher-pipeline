package <- snakemake@params[["selected_package"]]

is_package_available <- require(package)

if (!isTRUE(is_package_available)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    if (grepl(package, "T2T", fixed = TRUE)) {
        BiocManager::install("GenomeInfoDbData", update = FALSE)
        install.packages(package)
    } else {
        BiocManager::install(package, update = FALSE)
    }
    quit(save = "no")
}
