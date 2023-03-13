# package <- snakemake@params[["selected_package"]]
package <- "workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz"
print(grepl("BSgenome.T2T.CHM13.V2_1.0.0.tar.gz", package, fixed = TRUE, perl = FALSE))

is_package_available <- require(package)

if (!isTRUE(is_package_available)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    }
    if (grepl("BSgenome.T2T.CHM13.V2_1.0.0.tar.gz", package, fixed = TRUE, perl = FALSE)) {
        print("T2T")
        BiocManager::install("GenomeInfoDbData", update = FALSE)
        install.packages(package, repos = NULL, type = "source")
    } else {
        BiocManager::install(package, update = FALSE)
    }
    quit(save = "no")
}
