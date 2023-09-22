# package <- snakemake@params[["selected_package"]]
args <- commandArgs(TRUE)
package <- args[1]

# Check if the package is already available
is_package_available <- require(package, character.only = TRUE)

if (!isTRUE(is_package_available)) {
    # Ensure BiocManager is available since it will be needed regardless of the condition
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    }

    # Condition 1: Install custom tar.gz Bsgenome package named BSgenome.T2T.CHM13.V2_1.0.0.tar.gz
    if (grepl("BSgenome.T2T.CHM13.V2_1.0.0.tar.gz", package, fixed = TRUE, perl = FALSE)) {
        BiocManager::install("GenomeInfoDbData", update = FALSE)
        install.packages(package, repos = NULL, type = "source")

        # Condition 2: Install standard Bsgenome packages (hg38/hg19/mm10)
    } else if (package %in% c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Mmusculus.UCSC.mm10")) {
        BiocManager::install(package, update = FALSE)

        # Condition 3: Install a custom package using devtools
    } else {
        # Ensure devtools is installed
        if (!require("devtools", quietly = TRUE)) {
            install.packages("devtools", repos = "http://cran.us.r-project.org")
        }
        devtools::install(package, dependencies = TRUE, update = FALSE)
    }

    # Exit after installation, if desired
    # quit(save = "no")
}
