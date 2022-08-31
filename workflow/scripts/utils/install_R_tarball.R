package <- snakemake@input[["tarball"]]
is_package_available <- require(package)

if (!isTRUE(is_package_available)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("GenomeInfoDbData")
    install.packages(package)
    quit(save = "no")
}
