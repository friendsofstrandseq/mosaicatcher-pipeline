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

    # Set download options for better reliability
    options(timeout = 300)  # Increase timeout to 5 minutes
    options(download.file.method = "libcurl")

    # Retry function for robust installation
    install_with_retry <- function(install_func, max_attempts = 3) {
        for (attempt in 1:max_attempts) {
            tryCatch({
                install_func()
                cat("Installation successful on attempt", attempt, "\n")
                return(TRUE)
            }, error = function(e) {
                cat("Installation attempt", attempt, "failed:", conditionMessage(e), "\n")
                if (attempt < max_attempts) {
                    cat("Retrying in 5 seconds...\n")
                    Sys.sleep(5)
                } else {
                    cat("All attempts failed\n")
                    stop(e)
                }
            })
        }
        return(FALSE)
    }

    # Condition 1: Install custom tar.gz Bsgenome package named BSgenome.T2T.CHM13.V2_1.0.0.tar.gz
    if (grepl("BSgenome.T2T.CHM13.V2_1.0.0.tar.gz", package, fixed = TRUE, perl = FALSE)) {
        install_with_retry(function() {
            BiocManager::install("GenomeInfoDbData", update = FALSE)
            install.packages(package, repos = NULL, type = "source")
        })

        # Condition 2: Install standard Bsgenome packages (hg38/hg19/mm10)
    } else if (package %in% c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Mmusculus.UCSC.mm10")) {
        install_with_retry(function() {
            BiocManager::install(package, update = FALSE, ask = FALSE)
        })

        # Condition 3: Install a custom package using devtools
    } else {
        # Ensure devtools is installed
        if (!require("devtools", quietly = TRUE)) {
            install.packages("devtools", repos = "http://cran.us.r-project.org")
        }
        install_with_retry(function() {
            devtools::install(package, dependencies = TRUE, update = FALSE)
        })
    }

    # Exit after installation, if desired
    # quit(save = "no")
}
