source("workflow/scripts/plotting/plot-clustering.R")
plot.clustering(
    inputfile = ".tests/output_hg19/mosaiclassifier/sv_calls/RPE-BM510/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE.tsv",
    bin.bed.filename = "workflow/data/bin_200kb_all.bed",
    position.outputfile = "TEST_POSITION.pdf",
    chromosome.outputfile = "TEST_chr.pdf"
)