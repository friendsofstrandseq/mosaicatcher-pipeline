source("/g/korbel2/weber/workspace/mosaicatcher-update/workflow/scripts/plotting/plot-clustering.R")
traceback()
plot.clustering(
    inputfile = ".tests/data_CHR17_NEW/RPE-BM510/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv",
    bin.bed.filename = "workflow/data/bin_200kb_all.bed",
    position.outputfile = "TEST_POSITION.pdf",
    # chromosome.outputfile = "TEST_chr.pdf"
)