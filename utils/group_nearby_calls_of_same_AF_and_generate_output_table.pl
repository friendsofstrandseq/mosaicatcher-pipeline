#!/usr/bin/perl -w
use strict;


my $min_AF_for_grouping =0.03;
my $AF_simiarity_threshold = 0.6;
my $input_file = $ARGV[0];
my $safe_llr_to_ref = 30;
if (!$ARGV[0]) {
	print STDERR "inputfile missing\n";
	print STDERR "-> cluster into groups by chromosome, if similar AF (relative threshold=$AF_simiarity_threshold), if same SV event, if SVs are directly adjacent, and if at least one cell shared\n";
	print STDERR "Ignores SVs with VAF<$min_AF_for_grouping, which are never grouped\n";
	die;
}

#--- read input file
my (%MAIN_HAP, %PASS_FAIL, %Sv_call_haplotype, %Num_bins, %Group_ID, %STARTs, %ENDs, %SV_TYPEs, %CELLs, %AFs, %SEEN, %averaged_AF, %MAIN_SV_TYPE, %Llr_to_2nds, %Llr_to_refs);
print STDERR "Reading $input_file...\n";
open FH, "$input_file" or die;
while (<FH>) {
        chomp;
        next if ($_ =~ /^chrom/); #ignore first line
        my ($chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd, $af, $passfail) = split (/[\t ]+/, $_); #split by TAB/whitespace
	#print STDERR "$chrom, $start, $end, $sample, $cell, $strand_state_class, $scalar, $num_bins, $sv_call_name, $sv_call_haplotype, $sv_call_name_2nd, $sv_call_haplotype_2nd, $llr_to_ref, $llr_to_2nd, $af, $passfail\n";
	#next unless ($passfail =~ /PASS/); #ignore FAIL
	unless ($SEEN{$chrom}{$start}) {
		push (@{$STARTs{$chrom}}, $start);
        	push (@{$ENDs{$chrom}}, $end);
		$SEEN{$chrom}{$start}=1;
	}
        push (@{$PASS_FAIL{$chrom}{$start}}, $passfail); 
	push (@{$SV_TYPEs{$chrom}{$start}}, $sv_call_name);
	push (@{$CELLs{$chrom}{$start}}, $cell);
	push (@{$AFs{$chrom}{$start}}, $af);
	push (@{$Llr_to_refs{$chrom}{$start}}, $llr_to_ref);
	push (@{$Llr_to_2nds{$chrom}{$start}}, $llr_to_2nd);
	push (@{$Num_bins{$chrom}{$start}}, $num_bins);
	push (@{$Sv_call_haplotype{$chrom}{$start}}, $sv_call_haplotype);
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

#-- define majority haplotype
foreach my $chrom (sort keys %STARTs) {
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
                my (%count_abundance, $main_hap);
                foreach my $hapname (@{$Sv_call_haplotype{$chrom}{$STARTs{$chrom}[$i]}}) {
                        $count_abundance{$hapname}++;
                }
                ($main_hap) = sort {$count_abundance{$b}<=>$count_abundance{$a}} keys %count_abundance;
                $MAIN_HAP{$chrom}{$STARTs{$chrom}[$i]} = $main_hap;
        }
}


#--- cluster into groups by chromosome, if similar averaged AF, same primary SV type, directly adjacent call, and at least one cell shared
my $group_name = 0;
foreach my $chrom (sort keys %STARTs) {
	my $previous_shared = 0;
        for (my $i=0; $i<@{$STARTs{$chrom}}-1; $i++) {
		next if ($averaged_AF{$chrom}{$STARTs{$chrom}[$i]} < $min_AF_for_grouping); #ignore events with too low VAF
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
		print STDERR "Grouping/merging $chrom\t$STARTs{$chrom}[$i]-$ENDs{$chrom}[$i]|$STARTs{$chrom}[$i+1]-$ENDs{$chrom}[$i+1] (AF=$averaged_AF{$chrom}{$STARTs{$chrom}[$i]}|$averaged_AF{$chrom}{$STARTs{$chrom}[$i+1]}) (SV_type=$MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i]}) group_name=$group_name ...\n";	
		$previous_shared = 1;
		$Group_ID{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}=$group_name;
		$Group_ID{$chrom}{$STARTs{$chrom}[$i+1]}{$ENDs{$chrom}[$i+1]}=$group_name;
	}
}


#--------- Generating output
print STDERR "--\nGenerating Results_Output_Table\n";
print "chrom, start, end, num_bins, consensus_sv, consensus_sv_call_haplotype, llr_to_ref_max, llr_to_2nd_max, af, segments\n";

foreach my $chrom (sort keys %STARTs) {
        my $previous_shared = 0;
	my $last_group_ID;
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {
		if (exists ($Group_ID{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) {
			#next if ($seen{$Group_ID{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}});
			$last_group_ID = $Group_ID{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]};
			my (@myStarts, @myEnds, @myPassFail, $myLlr_to_refs);
			my $segments = "$chrom:$STARTs{$chrom}[$i]-$ENDs{$chrom}[$i]";
			push (@myStarts, $STARTs{$chrom}[$i]); 
			push (@myEnds, $ENDs{$chrom}[$i]);
			my $new_passfail;
			#print STDERR "tst\n";
			#print STDERR "---> $PASS_FAIL{$chrom}[$i][0]\n";
			#print STDERR "---> @{$PASS_FAIL{$chrom}[$i]}\n";
			($new_passfail) = sort {$b cmp $a} (@{$PASS_FAIL{$chrom}{$STARTs{$chrom}[$i]}});
			push (@myPassFail, $new_passfail);
			my $end =0;
			my $y=$i;
				#-- prepare consensus report for merged calls
			my ($af) = sort {$b <=> $a} @{$AFs{$chrom}{$STARTs{$chrom}[$i]}}; #use maximum in this consensus report
                        my ($llr_to_ref) = sort {$b <=> $a} @{$Llr_to_refs{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
                        my ($llr_to_2nd) = sort {$b <=> $a} @{$Llr_to_2nds{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
                        my ($num_bins) = sort {$b <=> $a} @{$Num_bins{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
			$myLlr_to_refs+= $llr_to_ref; #add up log likelihoods
			while (!$end) { #get other events falling onto same group ID
				$y++; #increment
				if ($y==@{$STARTs{$chrom}}) {
					$y--;
					$i=$y;
					last;
				}
				#print "($y)|", scalar (@{$STARTs{$chrom}}), " $chrom|$Group_ID{$chrom}{$STARTs{$chrom}[$y]}{$ENDs{$chrom}[$y]}|--->";
				if (exists ($Group_ID{$chrom}{$STARTs{$chrom}[$y]}{$ENDs{$chrom}[$y]})) {
					if ($last_group_ID ne $Group_ID{$chrom}{$STARTs{$chrom}[$y]}{$ENDs{$chrom}[$y]}) { #new group starts
						$y--;
                                        	$i=$y;
                                        	last;
					}
					push (@myStarts, $STARTs{$chrom}[$y]);                                         
                        		push (@myEnds, $ENDs{$chrom}[$y]);
					my $new_passfail;
                        		($new_passfail) = sort {$b cmp $a} @{$PASS_FAIL{$chrom}{$STARTs{$chrom}[$y]}};
                        		push (@myPassFail, $new_passfail);
					
					$segments .= "|$chrom:$STARTs{$chrom}[$y]-$ENDs{$chrom}[$y]";
					($af) = sort {$b <=> $a} ($af, @{$AFs{$chrom}{$STARTs{$chrom}[$i]}}); #use maximum in this consensus report
		                        ($llr_to_ref) = sort {$b <=> $a} ($llr_to_ref, @{$Llr_to_refs{$chrom}{$STARTs{$chrom}[$i]}}); #use max in this consensus report
               		                ($llr_to_2nd) = sort {$b <=> $a} ($llr_to_2nd, @{$Llr_to_2nds{$chrom}{$STARTs{$chrom}[$i]}}); #use max in this consensus report
                                        ($num_bins) = sort {$b <=> $a} ($num_bins, @{$Num_bins{$chrom}{$STARTs{$chrom}[$i]}}); #use max in this consensus report
					$myLlr_to_refs+= $llr_to_ref; #add up log likelihoods
				} else {
					$y--;
					$i=$y; #jump forward
					$end=1;
				}
			}
			#--> check here for PASS-FAIL in Groups using @myPassFail
			my $GroupPassFail;
			($GroupPassFail) = sort {$b cmp $a} @myPassFail;
			unless ($GroupPassFail =~ /PASS/) {
				next unless ($myLlr_to_refs >= $safe_llr_to_ref);
			}
			print "$chrom, $myStarts[0], $myEnds[-1], $num_bins, $MAIN_SV_TYPE{$chrom}{$myStarts[0]}, $MAIN_HAP{$chrom}{$myStarts[0]}, $myLlr_to_refs, N/A, $af, \[Group\_$Group_ID{$chrom}{$myStarts[0]}{$myEnds[0]}\/$segments\]\n";
			#print "$chrom, $myStarts[0], $myEnds[-1], $num_bins, $MAIN_SV_TYPE{$chrom}{$myStarts[0]]}\n";	
		} else {
			#-- prepare consensus report for singlish calls (calls not falling into a merge-group)
			my $passfail_single;
			($passfail_single) = sort {$a cmp $b} @{$PASS_FAIL{$chrom}{$STARTs{$chrom}[$i]}}; #remove singlish calls with a FAIL (introduced 31/10)
			next unless ($passfail_single =~ /PASS/); 
			my ($af) = sort {$b <=> $a} @{$AFs{$chrom}{$STARTs{$chrom}[$i]}}; #use maximum in this consensus report
			my ($llr_to_ref) = sort {$b <=> $a} @{$Llr_to_refs{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
			my ($llr_to_2nd) = sort {$b <=> $a} @{$Llr_to_2nds{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
			my ($num_bins) = sort {$b <=> $a} @{$Num_bins{$chrom}{$STARTs{$chrom}[$i]}}; #use max in this consensus report
			print "$chrom, $STARTs{$chrom}[$i], $ENDs{$chrom}[$i], $num_bins, $MAIN_SV_TYPE{$chrom}{$STARTs{$chrom}[$i]}, $MAIN_HAP{$chrom}{$STARTs{$chrom}[$i]}, $llr_to_ref, $llr_to_2nd, $af, [$chrom:$STARTs{$chrom}[$i]-$ENDs{$chrom}[$i]]\n";
			
		}
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

