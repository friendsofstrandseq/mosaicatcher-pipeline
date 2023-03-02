@files = <*.bam>;

foreach $file (@files) {
print "mv ".$file." ";
$file=~s/.sort.mdup.bam.sc_pre_mono_sort_for_mark_uniq.bam//;
print " ".$file.".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam\n\n";
print "samtools index ".$file.".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam\n\n";

}
