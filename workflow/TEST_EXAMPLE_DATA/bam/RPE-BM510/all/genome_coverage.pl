#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $path = $ARGV[0] or die;
my @bam_files = glob("$path*.bam");

open OUT, '>', 'genome_coverage.txt' or die "Can't read from file";

print OUT "Filename\tCoverage\tDepth_of_Coverage\n";

my %coverage = ();
my $length_map = 0;
my $pos_cov = 0;
my $bases = 0;

my $len = 2_745_186_691; #human_GRCh38 autosomes only

foreach my $bam (@bam_files) { 
warn "Working on $bam ...\n";
	foreach my $chr (1..22) {
		$chr = 'chr'.$chr;
		warn "Processing $chr ...\n";

		open F1, "samtools view -bh -q 10 -F1024 $bam $chr| samtools mpileup - |";

		while (<F1>) {
			chomp;
			my ($chr, $pos, $cov) = (split /\t/, $_)[0,1,3];
			next if $cov == 0;
			$bases += $cov;
			$pos_cov++;
			$coverage{$cov}++;
		}
	}

my $coverage_perc = ($pos_cov/$len) * 100;
my $depth_coverage = $bases/$len;

print OUT "$bam\t$coverage_perc\t$depth_coverage\n";

$pos_cov = 0;
$bases = 0;

}
