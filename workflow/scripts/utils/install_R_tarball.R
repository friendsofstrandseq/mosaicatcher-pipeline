is_snakemake_available <- require(snakemake)

if (isTRUE(is_snakemake_available)) {
    package <- snakemake@input[["tarball"]]
}   else {
    args <- commandArgs(trailingOnly = T)
    package <- args[1]
}

is_package_available <- require(package)

if (!isTRUE(is_package_available)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("GenomeInfoDbData")
    install.packages(package)
    quit(save = "no")
}
