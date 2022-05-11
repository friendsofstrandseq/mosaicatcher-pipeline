library(devtools)

load_all("/g/korbel2/weber/Gits/Rsamtools/")
load_all("/g/korbel2/weber/Gits/GenomicAlignments/")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")

bamfile <- "/g/korbel2/weber/MosaiCatcher_files/tmp/h/RPE1WTPE20401.sort.mdup.bam"
bamindex <- "/g/korbel2/weber/MosaiCatcher_files/tmp/h/RPE1WTPE20401.sort.mdup.bam.bai"
region <- vcf2ranges("/g/korbel2/weber/MosaiCatcher_output_sample_KG_chr21_test/snv_genotyping/RPE1-WT/chr21.vcf")

readGAlignmentPairs(bamfile, index = bamindex, param = Rsamtools::ScanBamParam(
    tag = "XA", which =
        range(region), what = c("seq", "qual", "mapq", "cigar", "flag"),
))