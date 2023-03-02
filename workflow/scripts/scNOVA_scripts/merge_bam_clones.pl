#!/usr/bin/perl


$input_subclonality=$ARGV[0];
$subclonality_colnames=$ARGV[1];
$line=$ARGV[2];
my $prefix = (split "scNOVA_input_user", $input_subclonality)[0];
chop($prefix);

if ( length($prefix) == 0 ){
    $prefix = ".";
}
print +($prefix)[0], "\n";
`mkdir -p $prefix/scNOVA_bam_merge`;

# `sed -n '1p' input_user/input_subclonality.txt > input_user/input_subclonality_colnames.txt`;
`sed -n '1p' $input_subclonality > $subclonality_colnames`;

@clone_all = ("clone1", "clone2", "clone3", "clone4", "clone5");
foreach $clone_ind (@clone_all) {

    # `grep "$clone_ind" input_user/input_subclonality.txt | wc -l > line.txt`;
    `grep "$clone_ind" $input_subclonality | wc -l > $line`;
    # open (FILE, "line.txt");
    open (FILE, $line);
	while (<FILE>){
	chomp;
	($x) = split("\t");
	print "$x\n";
	}
    close (FILE);

    if ($x > 0) {

		# `grep "$clone_ind" input_user/input_subclonality.txt > input_user/input_subclonality.$clone_ind.pre.txt`;
		`grep "$clone_ind" $input_subclonality > $prefix/scNOVA_input_user/input_subclonality.$clone_ind.pre.txt`;
		`cat $subclonality_colnames $prefix/scNOVA_input_user/input_subclonality.$clone_ind.pre.txt > $prefix/scNOVA_input_user/input_subclonality.$clone_ind.txt`;
		# `rm $prefix/scNOVA_input_user/input_subclonality.$clone_ind.pre.txt`; 
		`mkdir -p $prefix/scNOVA_bam_merge/$clone_ind`;

    	open (FILE, "$prefix/scNOVA_input_user/input_subclonality.$clone_ind.txt");
    	while (<FILE>) {
		chomp;
		($cell, $clone) = split("\t");

		print "cell: $cell\n";
		print "clone: $clone\n";
		print "---------\n";

		`cp $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam $prefix/scNOVA_bam_merge/$clone_ind`;
		`cp $prefix/scNOVA_bam_modified/$cell.sc_pre_mono_sort_for_mark_uniq.bam.bai $prefix/scNOVA_bam_merge/$clone_ind`;

		}
		close (FILE);
		`samtools merge -f $prefix/scNOVA_bam_merge/$clone_ind.merge.bam $prefix/scNOVA_bam_merge/$clone_ind/*.bam`;
		`samtools index $prefix/scNOVA_bam_merge/$clone_ind.merge.bam`;
		
    }

} 
