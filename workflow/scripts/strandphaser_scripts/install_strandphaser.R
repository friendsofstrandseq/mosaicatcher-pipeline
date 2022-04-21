# !/usr/bin/env Rscript


## #!/usr/bin/env Rscript
# Sys.setenv(Renv='PWD')

# Rscript version updated according to README 30-Mar-22 commit 69c9fb4
# Long time execution : ~20 min

# install.packages("devtools", repos = "http://cran.us.r-project.org")
library(devtools)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# BiocManager::install("GenomicRanges")
# BiocManager::install("GenomicAlignments")
# # withr::with_libpaths(new = "utils/R-packages", install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master"), "prefix")
# install_github("daewoooo/StrandPhaseR")


# Peter Ebert version

if (is.element("StrandPhaseR", installed.packages()[, 1])) {
    print("Removing previously installed version of StrandPhaseR")
    remove.packages("StrandPhaseR")
}

args <- commandArgs(trailingOnly = TRUE)

# Change PE to fixed commit
# git.commit <- "69c9fb4"
# git.repo <- "https://github.com/daewoooo/StrandPhaseR"
git.commit <- args[1]
git.repo <- args[2]

# Change from GIT to HTTPS
# if (is.na(git.repo)) {
#     git.repo <-
# }

devtools::install_git(
    git.repo,
    ref = git.commit,
    dependencies = FALSE,
    upgrade = FALSE
)

quit(save = "no")