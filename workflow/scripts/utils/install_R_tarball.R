package <- snakemake@input[["tarball"]]
is_package_available <- require(package)

if (!isTRUE(is_strandphaser_available)) {

    install.packages(package)
    quit(save = "no")
}
