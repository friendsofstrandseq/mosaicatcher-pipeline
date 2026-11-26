# package <- snakemake@params[["selected_package"]]
args <- commandArgs(TRUE)
package <- args[1]

# Check if the package is already available.
# `package` may be a tarball PATH (e.g. ".../BSgenome.Cfamiliaris.NCBI.canFam4_1.0.0.tar.gz") or a plain
# package name. Using the tarball path directly with require() never matches an installed
# package, so derive the real package name first (strip dir + "_<version>.tar.gz").
pkg_name <- if (grepl("\\.tar\\.gz$", package)) sub("_[0-9].*$", "", basename(package)) else package
is_package_available <- requireNamespace(pkg_name, quietly = TRUE)

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

    # Condition 1: Install a custom BSgenome source tarball (custom_bsgenome_tarball in
    # config, e.g. canFam4). Any local .tar.gz is installed from source.
    if (grepl("\\.tar\\.gz$", package)) {
        install_with_retry(function() {
            BiocManager::install("GenomeInfoDbData", update = FALSE)
            install.packages(package, repos = NULL, type = "source")
        })

        # Condition 2: Install standard BSgenome packages from Bioconductor
    # Matches: BSgenome.<Species>.UCSC.<assembly>
    # Examples: BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hs1 (T2T), BSgenome.Mmusculus.UCSC.mm39, etc.
    } else if (grepl("^BSgenome\\.[A-Z][a-z]+\\.UCSC\\.", package, perl = TRUE)) {
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
