#!/usr/bin/perl -w
use strict;

my $min_N_inv = 3;

#Use input file with format "chrom	start	end	sample	cell	class	scalar	num_bins	sv_call_name	sv_call_haplotype	sv_call_name_2nd	sv_call_haplotype_2nd	llr_to_ref	llr_to_2nd"
#
# e.g. 
# 

my $input_file = $ARGV[0];

if (!$ARGV[0]) {
	print STDERR "Input filname missinf=g (use e.g. sv_calls_txt_file_all/RPE1-WT/100000_fixed_norm.selected_j0.01_s0.1/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0.05_regfactor6.txt )\n";
	print STDERR "Filters inversions unless they are seen at least $min_N_inv times.\n";
	print STDERR "Filters deletions only seen in WW and CC chromosomes\n";
	print STDERR "Filters duplications only seen in WW and CC chromosomes (but gives inv-dups seen in such context a PASS)\n";
	exit;
}

print STDERR "Processing $input_file...\n";


#--- read input file
my (%STARTs, %ENDs, %SV_TYPEs, %STRAND_STATESs, %FILTER, %REPORTs);
open FH, "$input_file" or die;
while (<FH>) {
	chomp;
	next if ($_ =~ /^chrom/); #ignore first line
	my ($chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd) = split (/\t/, $_); #split by TAB
	push (@{$STARTs{$chrom}}, $start);
	push (@{$ENDs{$chrom}}, $end);
	push (@{$SV_TYPEs{$chrom}}, $sv_call_name);
	push (@{$STRAND_STATESs{$chrom}}, $strand_state_class);
	push (@{$REPORTs{$chrom}}, $_);
	}
close FH;


#--- 1. filter DEL calls only seen on WW or CC chromosomes
foreach my $chrom (sort keys %STARTs) {
	for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		next unless ($SV_TYPEs{$chrom}[$i] =~ /del_h[12]/); #only look at het-del calls
		$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'FAIL' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})); #initialize
		#print STDERR "$chrom $SV_TYPEs{$chrom}[$i]\n";
		$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($STRAND_STATESs{$chrom}[$i] eq 'WC' or $STRAND_STATESs{$chrom}[$i] eq 'CW');
		
	}
}

#--- 2. filter INV calls seen less than Min_N_inv times
foreach my $chrom (sort keys %STARTs) {
        my %inv_recurrence=();
	for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
                next unless ($SV_TYPEs{$chrom}[$i] =~ /inv_h[12]/); #only look at het-inv calls
                $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'FAIL' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})); #initialize
                #print STDERR "$chrom $SV_TYPEs{$chrom}[$i]\n";
		$inv_recurrence{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}++;
                $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($inv_recurrence{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $min_N_inv);

        }
}

#--- 3. filter DUP calls only seen on WW or CC chromosomes
foreach my $chrom (sort keys %STARTs) {
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
                next unless ($SV_TYPEs{$chrom}[$i] =~ /dup_h[12]/); #only look at het-dup calls
                $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'FAIL' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}) or $SV_TYPEs{$chrom}[$i] =~ /idup_h[12]/); #initialize
                #print STDERR "$chrom $SV_TYPEs{$chrom}[$i]\n";
                $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($STRAND_STATESs{$chrom}[$i] eq 'WC' or $STRAND_STATESs{$chrom}[$i] eq 'CW');

        }
}


#--- output file with filters
print "chrom\tstart\tend\tsample\tcell\tstrand_state_class\tscalar\tnum_bins\tsv_call_name\tsv_call_haplotype\tsv_call_name_2nd\tsv_call_haplotype_2nd\tllr_to_ref\tllr_to_2nd\tpass/fail\n";
foreach my $chrom (sort keys %STARTs) {
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'PASS' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})); #initialize
		print "$REPORTs{$chrom}[$i]\t", $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}, "\n";
	}
}
