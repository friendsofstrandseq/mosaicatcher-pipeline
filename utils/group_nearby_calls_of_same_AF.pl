#!/usr/bin/perl -w
use strict;


my $AF_simiarity_threshold = 0.25;
my $input_file = $ARGV[0];

if (!$ARGV[0]) {
	print STDERR "inputfile missing\n";
	print STDERR "-> cluster into groups by chromosome, if similar AF (relative threshold=$AF_simiarity_threshold), if same SV event, if SVs are directly adjacent, and if at least one cell shared\n";
	die;
}

#--- read input file
my (%STARTs, %ENDs, %SV_TYPEs, %CELLs, %AFs, %SEEN, %averaged_AF, %MAIN_SV_TYPE);
print STDERR "Reading $input_file...\n";
open FH, "$input_file" or die;
while (<FH>) {
        chomp;
        next if ($_ =~ /^chrom/); #ignore first line
        my ($chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd, $af, $passfail) = split (/[\t ]+/, $_); #split by TAB/whitespace
	#print STDERR "$chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd, $af, $passfail\n";
	next unless ($passfail =~ /PASS/); #ignore FAIL
	unless ($SEEN{$chrom}{$start}) {
		push (@{$STARTs{$chrom}}, $start);
        	push (@{$ENDs{$chrom}}, $end);
		$SEEN{$chrom}{$start}=1;
	}
        push (@{$SV_TYPEs{$chrom}{$start}}, $sv_call_name);
	push (@{$CELLs{$chrom}{$start}}, $cell);
	push (@{$AFs{$chrom}{$start}}, $af);
}
%SEEN=(); #reinitialize
close FH;


#-- compute average AF per event
foreach my $chrom (sort keys %STARTs) {
	#print STDERR "test\n";
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		$averaged_AF{$chrom}{$STARTs{$chrom}[$i]} = &average (@{$AFs{$chrom}{$STARTs{$chrom}[$i]}});	
	}
}

#-- define majority SV type (haplotype-independent)
foreach my $chrom (sort keys %STARTs) {
	for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		my (%count_abundance, $main_SV_type);
		foreach my $callname (@{$SV_TYPEs{$chrom}{$STARTs{$chrom}[$i]}}) {
			$callname =~ s/_h[12]//;
			#print STDERR "$callname\n";
			$count_abundance{$callname}++;
		}
		($main_SV_type) = sort {$count_abundance{$b}<=>$count_abundance{$a}} keys %count_abundance;
		$MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i]} = $main_SV_type;	
	}
}

#--- cluster into groups by chromosome, if similar averaged AF, same primary SV type, directly adjacent call, and at least one cell shared
my $group_name = 0;
foreach my $chrom (sort keys %STARTs) {
	my $previous_shared = 0;
        for (my $i=0; $i<@{$STARTs{$chrom}}-1; $i++) {
		unless ($ENDs{$chrom}[$i] == $STARTs{$chrom}[$i+1]) { #only consider events that are directly adjacent
			$previous_shared=0;#initialize
			next;
		}
		unless (($averaged_AF{$chrom}{$STARTs{$chrom}[$i]}/$averaged_AF{$chrom}{$STARTs{$chrom}[$i+1]}>= (1-$AF_simiarity_threshold)) and ($averaged_AF{$chrom}{$STARTs{$chrom}[$i+1]}/$averaged_AF{$chrom}{$STARTs{$chrom}[$i]}>= (1-$AF_simiarity_threshold))) { #compare AFs
			$previous_shared=0;#initialize
			next;
		}
		unless ($MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i]} eq $MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i+1]}) { #must be same primary SV type
			$previous_shared=0;#initialize
			next;
		}
		my $shared_cells=0;
		foreach my $cell1 (@{$CELLs{$chrom}{$STARTs{$chrom}[$i]}}) {
			last if ($shared_cells==1);
			foreach my $cell2 (@{$CELLs{$chrom}{$STARTs{$chrom}[$i+1]}}) {
				$shared_cells = 1 if ($cell1 eq $cell2);
			}
		}
		if ($shared_cells == 0) {
			$previous_shared=0;#initialize
			next;
		}
		$group_name++ unless ($previous_shared);
		print "SHARED: $chrom\t$STARTs{$chrom}[$i]-$ENDs{$chrom}[$i]|$STARTs{$chrom}[$i+1]-$ENDs{$chrom}[$i+1] (AF=$averaged_AF{$chrom}{$STARTs{$chrom}[$i]}|$averaged_AF{$chrom}{$STARTs{$chrom}[$i+1]}) (SV_type=$MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i]}) group_name=$group_name\n";	
		$previous_shared = 1;
	}
}

#---------------
#--- subroutines
#---------------
sub average {
	my $av=0;
	foreach my $cnt (@_) {
		$av+=$cnt/scalar(@_);
	}
	return $av;
}

