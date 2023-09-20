package <- snakemake@params[["selected_package"]]

is_package_available <- require(package, character.only = TRUE)

if (!isTRUE(is_package_available)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    }

    if (grepl("BSgenome.T2T.CHM13.V2_1.0.0.tar.gz", package, fixed = TRUE, perl = FALSE)) {
        print("T2T")
        BiocManager::install("GenomeInfoDbData", update = FALSE)
        install.packages(package, repos = NULL, type = "source")
    } else {
        # Check if devtools is installed, if not install it
        if (!require("devtools", quietly = TRUE)) {
            install.packages("devtools", repos = "http://cran.us.r-project.org")
        }
        # Use devtools::install to install the package
        devtools::install(package, dependencies = TRUE, update = FALSE)
    }
    quit(save = "no")
}
