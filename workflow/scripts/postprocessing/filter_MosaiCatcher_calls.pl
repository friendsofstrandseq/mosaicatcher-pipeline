#!/usr/bin/perl -w
use strict;

my $min_N_inv = 3;
my $min_WC = 1/3;
my $safe_llr_to_ref = 30;
# my $SegDup_file = "../../data/segdups/segDups_hg38_UCSCtrack.bed.gz";
my $SegDup_file = $ARGV[1];
my $MaxSegDup_overlap = 0.5;

#Use input file with format "chrom	start	end	sample	cell	class	scalar	num_bins	sv_call_name	sv_call_haplotype	sv_call_name_2nd	sv_call_haplotype_2nd	llr_to_ref	llr_to_2nd	af"
#
# e.g. 
# 

my $input_file = $ARGV[0];

if (!$ARGV[0]) {
	print STDERR "Input filname missinf=g (use e.g. sv_calls_txt_file_all/RPE1-WT/100000_fixed_norm.selected_j0.01_s0.1/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0.05_regfactor6.txt )\n";
	print STDERR "Filters inversions unless they are seen at least $min_N_inv times.\n";
	print STDERR "Filters deletions seen in not more than $min_WC WC chromosomes\n";
	print STDERR "Filters duplications seen in not more than $min_WC WC chromosomes (but gives inv-dups seen in such context a PASS)\n";
	print STDERR "Del and Dup events with llr_to_ref>=$safe_llr_to_ref will never be masked\n";
	printf STDERR "uses SegDup file $SegDup_file and removes all Dels overlapping with SegDups by >%4.1f percent \n", $MaxSegDup_overlap*100;
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

#-- 4. read SegDup table
print STDERR "Parsing $SegDup_file... ";
open FH_1, "gunzip -c $SegDup_file |" or die;
my (%SegDupStart, %SegDupEnd);
my $cntS;
my %seenSegDup;
while (<FH_1>) {
	chomp;
	next if ($_=~/^#/); #ignore first line
	my (undef, $chr, $start, $end) = split (/[\t ]+/, $_);
	next if (exists ($seenSegDup{$chr}{$start}{$end}));
	$seenSegDup{$chr}{$start}{$end}=1;
	#print STDERR "($chr, $start, $end)\n";
	push (@{$SegDupStart{$chr}}, $start);
	push (@{$SegDupEnd{$chr}}, $end);
	$cntS++;
		
}
print STDERR "$cntS events parsed.\n";
close FH_1;

#-- 5. filter Del events overlapping SegDups to much (based on defined value)
print STDERR "Testing for overlap with SegDups...\n";
#my %SEGDUPOV;
my $UCSC_artefact_chr="chr2";
my $UCSC_artefact_start=90400000;
my $UCSC_artefact_end=91400000;

foreach my $chrom (sort keys %STARTs) {
	my %Seen;
        for (my $i=0; $i<@{$STARTs{$chrom}}; $i++) {	
		next unless ($SV_TYPEs{$chrom}[$i] =~ /del/); #only look at del calls
		my $stD=$STARTs{$chrom}[$i]; my $enD=$ENDs{$chrom}[$i];
		next if (exists($Seen{$stD}{$enD}));
		$Seen{$stD}{$enD}=1;
		my $overlap_string=sprintf "0"x($enD-$stD+1); #Initialize" make long string of "0" with size ($enD-$stD+1)	
		for (my $y=0; $y<@{$SegDupStart{$chrom}}; $y++) {
			my $stS=$SegDupStart{$chrom}[$y]; my $enS=$SegDupEnd{$chrom}[$y];
			if (($stD<=$stS) && ($stS<=$enD)) {
				my ($EndPos) = sort {$a<=>$b} ($enD, $enS); #sort smaller EndVal
				my $DupOvL = $EndPos-$stS+1;
				#print "$overlap_string\n";
				#print "vorher:", length ($overlap_string), "\n";
				substr ($overlap_string, $stS-$stD+1, $DupOvL) = sprintf "1"x$DupOvL; #fill in "1" for each overlapping base	
			} elsif (($stS<=$stD) && ($stD<=$enS)) {
				my ($EndPos) = sort {$a<=>$b} ($enD, $enS); #sort smaller EndVal
				my $DupOvL = $EndPos-$stD+1;
				#print "vorher:", length ($overlap_string), "\n";
				substr ($overlap_string, 0, $DupOvL) = sprintf "1"x$DupOvL; #fill in "1" for each overlapping base 
				#print "nachher:", length ($overlap_string), "\n";
				
			}
			
		}
		my $overlap=0;
                for (my $M=0; $M<length($overlap_string); $M++) { 
                     $overlap+=1 if (substr($overlap_string, $M, 1) eq "1");
                }
		$overlap/=length($overlap_string);
		#$SEGDUPOV{$chrom}{$stD}{$enD}=$overlap;
		$FILTER{$chrom}{$stD}{$enD}= sprintf ("FAIL(SegDup:%4.2f)", $overlap) if ($overlap > $MaxSegDup_overlap);		
		#------------------------------------------------------------
		#-- remove build38 / ucsc 'hole' region in annotated assembly
		#------------------------------------------------------------
		next unless ($chrom eq $UCSC_artefact_chr);
		my $UCSC_exception=0;
		if (($stD<=$UCSC_artefact_start) && ($UCSC_artefact_start<=$enD)) {
                                my ($EndPos) = sort {$a<=>$b} ($enD, $UCSC_artefact_end); #sort smaller EndVal
                                my $DupOvL = $EndPos-$UCSC_artefact_start+1;
                                substr ($overlap_string, $UCSC_artefact_start-$stD+1, $DupOvL) = sprintf "1"x$DupOvL; #fill in "1" for each overlapping base
                		$UCSC_exception=1; 
		} elsif (($UCSC_artefact_start<=$stD) && ($stD<=$UCSC_artefact_end)) {
                       my ($EndPos) = sort {$a<=>$b} ($enD, $UCSC_artefact_end); #sort smaller EndVal
                       my $DupOvL = $EndPos-$stD+1;
                       substr ($overlap_string, 0, $DupOvL) = sprintf "1"x$DupOvL; #fill in "1" for each overlapping base
                	$UCSC_exception = 1;
		}
		next unless ($UCSC_exception);
		for (my $M=0; $M<length($overlap_string); $M++) {
                     $overlap+=1 if (substr($overlap_string, $M, 1) eq "1");
                }
		$overlap/=length($overlap_string);
		$FILTER{$chrom}{$stD}{$enD}= sprintf ("FAIL(SegD/UCSC:%4.2f)", $overlap) if ($overlap > $MaxSegDup_overlap);
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
			
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= sprintf ("FAIL(WC:%4.2f)", $iterate) if ($iterate <= $min_WC);
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= "PASS" if ($iterate > $min_WC or $LLR_TO_REF{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $safe_llr_to_ref);
		   }
		   if (exists ($TESTED_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]})) { #check DUP to identify events seen not largely in CC and WW chromosomes
                        my $iterate=0;
                        foreach my $val (@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}}) {
                                $iterate+=$val/@{$FILTER_ARR_DUP{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}};
                        }
                        $FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= sprintf ("FAIL(WC:%4.2f)", $iterate) if ($iterate <= $min_WC);
			$FILTER{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]}= "PASS" if ($iterate > $min_WC or $LLR_TO_REF{$chrom}{$STARTs{$chrom}[$i]}{$ENDs{$chrom}[$i]} >= $safe_llr_to_ref);
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

#--- subroutines
#sub compare_overlap {
#	my ($chr, $stD, $enD, $stS, $enS) = @_;
#	my $flag=0;
#	my $SV_length = $enD-$stD+1;
#	my $SegDup_length = $enS-$stS+1;
#	my $overlap_length=0;
#	my $Rel_overlap=0;
#	#-- assess_overlap
#	if (($enS>=$stD) and ($enS<=$enD) and ($stS>=$stD)) { #segDup is fully contained in event
#		$overlap_length = $enS-$stS+1;
#		$Rel_overlap = $overlap_length/$SV_length;
#		print  "($chr: $stD, $enD, $stS, $enS) -> overlap: $overlap_length (SV_length=$SV_length; SegDup_length=$SegDup_length) | relative_overlap=$Rel_overlap\n" if ($Rel_overlap>=$MaxSegDup_overlap);
#		$RelSegDupOverlap{$chr}{$stD}{$enD} = 1 if ($Rel_overlap>=$MaxSegDup_overlap);	
#	}
#	$flag = 1 if ($Rel_overlap>=$MaxSegDup_overlap);
#	print "--->$flag\n" if ($flag);
#	return $flag;
#}
