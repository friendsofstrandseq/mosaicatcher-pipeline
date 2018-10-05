#!/usr/bin/perl -w
use strict;

my $min_N_inv = 3;
my $min_WC = 0.2;
my $safe_llr_to_ref = 20;

#Use input file with format "chrom	start	end	sample	cell	class	scalar	num_bins	sv_call_name	sv_call_haplotype	sv_call_name_2nd	sv_call_haplotype_2nd	llr_to_ref	llr_to_2nd	af"
#
# e.g. 
# 

my $input_file = $ARGV[0];

if (!$ARGV[0]) {
	print STDERR "Input filname missinf=g (use e.g. sv_calls_txt_file_all/RPE1-WT/100000_fixed_norm.selected_j0.01_s0.1/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0.05_regfactor6.txt )\n";
	print STDERR "Filters inversions unless they are seen at least $min_N_inv times.\n";
	print STDERR "Filters deletions seen in less then $min_WC WC chromosomes\n";
	print STDERR "Filters duplications seen in less then $min_WC WC chromosomes (but gives inv-dups seen in such context a PASS)\n";
	print STDERR "Del and Dup events with llr_to_ref>=$safe_llr_to_ref will never be masked\n";
	exit;
}

print STDERR "Processing $input_file...\n";


#--- read input file
my (%STARTs, %ENDs, %SV_TYPEs, %STRAND_STATESs, %FILTER, %REPORTs, %FILTER_ARR, %FILTER_ARR_DUP, %TESTED_DEL, %TESTED_DUP, %inv_recurrence, %LLR_TO_REF);
open FH, "$input_file" or die;
while (<FH>) {
	chomp;
	next if ($_ =~ /^chrom/); #ignore first line
	my ($chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd, $af) = split (/\t/, $_); #split by TAB
	push (@{$STARTs{$chrom}}, $start);
	push (@{$ENDs{$chrom}}, $end);
	push (@{$SV_TYPEs{$chrom}}, $sv_call_name);
	push (@{$STRAND_STATESs{$chrom}}, $strand_state_class);
	#$_ =~ /[\t ]*/\t/;
	push (@{$REPORTs{$chrom}}, $_);

	# Del and Dup with high llr will never be filtered. Record highest llr for event at certain position
	$LLR_TO_REF{$chrom}{$start}{$end} = $llr_to_ref unless (exists ($LLR_TO_REF{$chrom}{$start}{$end}));
	$llr_to_ref = 100 if ($llr_to_ref =~ /Inf/);
	$LLR_TO_REF{$chrom}{$start}{$end} = $llr_to_ref if ($llr_to_ref > $LLR_TO_REF{$chrom}{$start}{$end});
	}
close FH;


#--- 1. filter DEL calls seen primarily on WW or CC chromosomes
foreach my $chrom (sort keys %STARTs) {
	for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		next unless ($SV_TYPEs{$chrom}[$i] =~ /del_h[12]/); #only look at het-del calls
		#FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'FAIL' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})); #initialize
		#print STDERR "$chrom $SV_TYPEs{$chrom}[$i]\n";
		if ($STRAND_STATESs{$chrom}[$i] eq 'WC' or $STRAND_STATESs{$chrom}[$i] eq 'CW') {
			push (@{$FILTER_ARR{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}, 1);
		} else {
			push (@{$FILTER_ARR{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}, 0);
		}
		$TESTED_DEL{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}=1;
		
	}
}

#--- 2. filter INV calls seen less than Min_N_inv times
foreach my $chrom (sort keys %STARTs) {
	for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
                next unless ($SV_TYPEs{$chrom}[$i] =~ /inv_h[12o]/); #only look at het-inv or hom-inv calls
                #$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'FAIL' unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})); #initialize
                #print STDERR "$chrom $SV_TYPEs{$chrom}[$i]\n";
		$inv_recurrence{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}++;
                #$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($inv_recurrence{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $min_N_inv);

        }
}

#--- 3. filter DUP calls seen primarily on WW or CC chromosomes
foreach my $chrom (sort keys %STARTs) {
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
                $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($SV_TYPEs{$chrom}[$i] =~ /^idup_h[12]/); #inverted duplications trigger 'pass'
                next unless ($SV_TYPEs{$chrom}[$i] =~ /^dup_h[12]/); #look at het-dup calls
                #$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = 'PASS' if ($STRAND_STATESs{$chrom}[$i] eq 'WC' or $STRAND_STATESs{$chrom}[$i] eq 'CW');
		
		if ($STRAND_STATESs{$chrom}[$i] eq 'WC' or $STRAND_STATESs{$chrom}[$i] eq 'CW') {
                        push (@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}, 1);
                } else {
                        push (@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}, 0);
                }
                $TESTED_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}=1;

        }
}


#--- output file with filters
print "chrom\tstart\tend\tsample\tcell\tstrand_state_class\tscalar\tnum_bins\tsv_call_name\tsv_call_haplotype\tsv_call_name_2nd\tsv_call_haplotype_2nd\tllr_to_ref\tllr_to_2nd\taf\tpass/fail\n";
foreach my $chrom (sort keys %STARTs) {
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		if (!exists ($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) {
		

		  if (exists ($TESTED_DEL{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) { #check DELs to identify events seen not largely in CC and WW chromosomes
			my $iterate=0;
			foreach my $val (@{$FILTER_ARR{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}) {
				$iterate+=$val/@{$FILTER_ARR{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}};
			}
			
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= sprintf ("FAIL(WC:%4.2f)", $iterate) if ($iterate < $min_WC);
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= "PASS" if ($iterate >= $min_WC or $LLR_TO_REF{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $safe_llr_to_ref);
		   }
		   if (exists ($TESTED_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) { #check DUP to identify events seen not largely in CC and WW chromosomes
                        my $iterate=0;
                        foreach my $val (@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}) {
                                $iterate+=$val/@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}};
                        }
                        $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= sprintf ("FAIL(WC:%4.2f)", $iterate) if ($iterate < $min_WC);
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= "PASS" if ($iterate >= $min_WC or $LLR_TO_REF{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $safe_llr_to_ref);
                   }
		}		
		unless (exists($FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) {
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= 'PASS';
			if (exists ($inv_recurrence{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) {
				$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} = sprintf ("FAIL(#inv:" . $inv_recurrence{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} . ")") if ($inv_recurrence{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} < $min_N_inv);	
			}	
		}
		print "$REPORTs{$chrom}[$i]\t", $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}, "\n";
	}
}
