@files = <*uniq.bam>;

print "`multiBamSummary BED-file --BED /g/korbel2/StrandSeq/Test_HJ/reference/bin_Genes_for_CNN_sort.txt --bamfiles ";
foreach $file (@files)
{
print  "$file ";
}
print "--extendReads --outRawCounts Deeptool_Genes_for_CNN_HG00512.tab -out Deeptool_Genes_for_CNN_HG00512.npz`;\n\n";

