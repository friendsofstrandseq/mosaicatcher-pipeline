#!/usr/bin/perl

$strandphaser_output=$ARGV[0];
$nucleosome_sampleA=$ARGV[1];
$nucleosome_sampleB=$ARGV[2];
$strandphaser_output_copy=$ARGV[3];

# WARNING : change input_user for strandphaser
my $prefix = (split "strandphaser", $strandphaser_output)[0];
chop($prefix);

if ( length($prefix) == 0 ){
    $prefix = ".";
}

# `sed -n '1p' input_user/strandphaser_output.txt > input_user/strandphaser_colnames.txt`;
`sed -n '1p' $strandphaser_output > $prefix/scNOVA_input_user/strandphaser_colnames.txt`;
`mkdir -p $prefix/nucleosome_sampleA`;
`mkdir -p $prefix/nucleosome_sampleB`;

@chrom_all = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX");

foreach $chrom_ind (@chrom_all) {

		# `grep "$chrom_ind\t" input_user/strandphaser_output.txt > input_user/strandphaser_output.$chrom_ind.pre.txt`;
		`grep "$chrom_ind\t" $strandphaser_output > $prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.pre.txt`;
		`cat $prefix/scNOVA_input_user/strandphaser_colnames.txt $prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.pre.txt > $prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.txt`;
		# `rm $prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.pre.txt`;

		print "chrom: $chrom_ind\n";


		# WARNING replace bam by scNOVA_bam
		`mkdir -p $prefix/scNOVA_bam_modified/$chrom_ind.H1`;
		`mkdir -p $prefix/scNOVA_bam_modified/$chrom_ind.H2`;
		print "$prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.txt\n";
		open (FILE, "$prefix/scNOVA_input_user/strandphaser_output.$chrom_ind.txt");
		while (<FILE>) {
		chomp;
		($sample, $cell, $chrom, $start, $end, $class, $hap1.cis.simil, $hap1.trans.simil, $hap2.cis.simil, $hap2.trans.simil) = split("\t");
  	$cell=~s/.sort.mdup.bam//;
		if ($chrom eq $chrom_ind){

			print "cell: $cell\n";
			print "chrom: $chrom\n";
			print "start: $start\n";
			print "end: $end\n";
			print "class: $class\n";
			print "---------\n";

			`samtools view -b $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.bam "$chrom:$start-$end" > $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam`;
			`samtools view -b $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.bam "$chrom:$start-$end" > $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam`;

			if ($class eq "WC") {
				`mv $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam $prefix/scNOVA_bam_modified/$chrom_ind.H1/`;
				`mv $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam $prefix/scNOVA_bam_modified/$chrom_ind.H2/`;
			} 

			if ($class eq "CW") {
				`mv $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam $prefix/scNOVA_bam_modified/$chrom_ind.H1/`;
				`mv $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam $prefix/scNOVA_bam_modified/$chrom_ind.H2/`;
			}

		}
		}	
		close (FILE);

	`find $prefix/scNOVA_bam_modified/$chrom_ind.H1/ -size 0 -print -delete`;
	`samtools merge $prefix/scNOVA_bam_modified/$chrom_ind.H1.bam $prefix/scNOVA_bam_modified/$chrom_ind.H1/*.seg.bam`;
	`samtools index $prefix/scNOVA_bam_modified/$chrom_ind.H1.bam`;
	# `rm -r $prefix/scNOVA_bam_modified/$chrom_ind.H1`;

	`find $prefix/scNOVA_bam_modified/$chrom_ind.H2/ -size 0 -print -delete`;
	`samtools merge $prefix/scNOVA_bam_modified/$chrom_ind.H2.bam $prefix/scNOVA_bam_modified/$chrom_ind.H2/*.seg.bam`;
	`samtools index $prefix/scNOVA_bam_modified/$chrom_ind.H2.bam`;
	# `rm -r $prefix/scNOVA_bam_modified/$chrom_ind.H2`;

	close (Chr_FILE);

	}

	# `samtools merge nucleosome_sampleA/result.H1.bam bam/*.H1.bam`;
	`samtools merge $nucleosome_sampleA $prefix/scNOVA_bam_modified/*.H1.bam`;
	# `samtools merge nucleosome_sampleB/result.H2.bam bam/*.H2.bam`;
	`samtools merge $nucleosome_sampleB $prefix/scNOVA_bam_modified/*.H2.bam`;
	# `samtools index nucleosome_sampleA/result.H1.bam`;
	`samtools index $nucleosome_sampleA`;
	# `samtools index nucleosome_sampleB/result.H2.bam`;
	`samtools index $nucleosome_sampleB`;
	# `cp input_user/strandphaser_output.txt input_user/strandphaser_output_copy.txt`;
	`cp $strandphaser_output $strandphaser_output_copy`;
