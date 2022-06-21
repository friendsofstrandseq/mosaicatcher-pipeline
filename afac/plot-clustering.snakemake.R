source("utils/plot-clustering.R")
plot.clustering(
    inputfile = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/sv_calls/RPE1-WT/100000.selected_j0.1_s0.5_scedist20/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt",
    bin.bed.filename = "utils/bin_200kb_all.bed",
    position.outputfile = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/sv_calls/RPE-BM510/100000.selected_j0.1_s0.5_scedist20/plots/sv_clustering/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE-position.pdf",
    chromosome.outputfile = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/sv_calls/RPE-BM510/100000.selected_j0.1_s0.5_scedist20/plots/sv_clustering/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE-chromosome.pdf"
)