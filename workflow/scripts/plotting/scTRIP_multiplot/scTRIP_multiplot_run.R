library(scTRIPmultiplot)
args <- commandArgs(TRUE)


counts_path <- args[1]
haplo_path <- args[2]
sv_path <- args[3]
chromosome <- args[4]
cell_id <- args[5]
savepath <- args[6]


scTRIPmultiplot::generate_multiplot(
    counts_path = counts_path,
    haplo_path = haplo_path,
    chromosome = chromosome,
    cell_id = cell_id,
    savepath = savepath,
    sv_path = sv_path,
    size = 1
)
