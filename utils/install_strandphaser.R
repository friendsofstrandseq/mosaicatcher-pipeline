#!/usr/bin/env Rscript
Sys.setenv(Renv='PWD')
library(devtools)
withr::with_libpaths(new = "utils/R-packages", install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master"), "prefix")
