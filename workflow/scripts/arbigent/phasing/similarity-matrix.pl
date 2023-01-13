#!/usr/bin/perl
#|
#|	similarity-matrix.pl
#|
#|	written by Connor Burbridge
#|
#|	This script has been written to compare samples in a vcf file. In order to find 
#|	the differences and similarities this script will look through the vcf and 
#|	calculate the ratios of identical calls between all samples. It will then create
#|	a tab-delimited matrix that allows for easy visualization of these ratios. There
#|	is an option for this matrix to have either the calculated percentage of these
#|	ratios or the raw counts of both samples. My original testing comments have been
#|  left in place for future testing purposes.
#|
#|	example call:
#|
#|  ./similarity-matrix.pl -c on -f ./mono-allelic-LR-26.vcf -o ./matrix.txt
#|
#|  Necessary arguments (in order):
#|
#|	-c -> either "on" or "off" 
#|		-the off option will leave the matrix with the raw identical call counts
#|		-the on option will calculate the ratios in percentages for the matrix
#|	-f -> the path to the input vcf file you want to use
#|	-o -> the path and or name of the matrix created by the script
#|__________________________________________________________________________________
####################################################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
my %options=();
getopts("c:f:o:", \%options);

my $vcffile = $options{f} if defined $options{f};
my $outfile = $options{o} if defined $options{o};
my $flag = $options{c} if defined $options{c};

unless (defined $vcffile)
{
	die "ERROR: vcf file not supplied! Please try again with an appropriate file. Documentation can be viewed within this script.\n";
}

unless ((defined $outfile) and !($outfile =~ /\s/))
{
	die "ERROR: outfile not supplied or contains white space! Please try again with a new outfile! Documentation can be viewed within this script.\n";
}

unless ($flag =~ /on|off/)
{
	die "ERROR: invalid option for -c was given! Please try again with \"on\" or \"off\"\nDocumentation can be viewed within this script.\n";
}

open INFILE, "$vcffile" or die "ERROR: could not open vcf file!\n";

my $line;
my $sampleindex = 0;
my @samples;
my @pass;
my %variants;
my $markername;
my $genotype;
my @header;
my @call;
my @marker;
my @splitsample;
my $samplecount;
my $sample;
my $current;
my $filelength = `grep -cv '^#' $vcffile`;
my $percent;
my $subpercent;
my $i;
#print $samplecount, "\n";
#print @samples, "\n";

while ($line = <INFILE>) #and ($line !~ /LcChr1\t3483114/))
{
	chomp $line;
	$percent = ($sampleindex/$filelength)*100;
	$subpercent = substr($percent,0,index($percent,'.') + 1 + 2);
	#print "\rBuilding hash: ", $subpercent, "%";
	#print $file, "\n";

	if ($line =~ /#CHROM/)
	{
		chomp $line;
		@header = split(/\t/, $line);
		for ($i = 9; $i < scalar(@header); $i = $i + 1)
		{
			push(@samples, $header[$i]);
			#print $header[$i], "\n";
		}
		$samplecount = @samples;
	}

	if ($line !~ /^#/)
	{
		chomp $line;
		@marker = split(/\t/, $line);
		$markername = $marker[0] . "p" . $marker[1];
		$sampleindex = $sampleindex + 1;
		#print $position, "\t", $genotype, "\n";

		foreach $sample (@samples)
		{
			#$variants{$sample} = {};
			#print "sample ", $sample, "\n";
            #print "samples ", @samples, "\n";
	
			for ($i = 8; $i <= ($samplecount + 8); $i = $i + 1)
			{
				@splitsample = $header[($i)]; #split('+', $header[($i)]);
				$current = $header[($i)];
				#print $current, "\n";
				#print "samplecount ", $samplecount, "\n";
				#print "sample and current ", $sample, ' ' , $current, "\n";
				#print "testres: ", $sample eq $current, "\n";
				if ($sample eq $current)
				{
					@call = split(/:/, $marker[$i]);
					$genotype = $call[0];
					#print $sample, "\t", $marker[$i], "\t", $markername, "\t", $genotype, "\n";
					#print "sup\n";
					$variants{$sample}{$markername} = $genotype;

				}
			}
		}
	}
}

print "\nHash table is ready. Comparing calls now...\n";
#print Dumper(%variants);
close INFILE or die "ERROR: could not close vcf file!\n";
open INFILE, $vcffile or die "ERROR: could not open vcf file!\n";
open OUTFILE, ">$outfile" or die "ERROR: could not open outfile!\n";
my $key;
my @keyline;
my $matchcount = 0;
my @vcfcalls;
my @hashcalls;
my $misscount = 0;
my $call;
my $totalcount = 0;
my $marker;
$sampleindex = 0;

print OUTFILE "sample_name";

foreach $sample (sort keys %variants)
{
	print OUTFILE "\t$sample";
}
print OUTFILE "\n";


foreach my $sample1 (sort keys %variants)
{
	$sampleindex = $sampleindex + 1;
	$percent = ($sampleindex/$samplecount)*100;
	$subpercent = substr($percent,0,index($percent,'.') + 1 + 2);
	print "\rComparing sample ", $sampleindex, " of ", $samplecount, ": ", $subpercent, "% done...";
	print OUTFILE $sample1;

	foreach my $sample2 (sort keys %variants)
	{
		$matchcount = 0;
		$totalcount = 0;
		foreach $marker (sort keys %{$variants{$sample2}})
		{
			#$call = $variants{$sample2}{$marker};
			if (($variants{$sample1}{$marker} =~ /\.\/\./) or ($variants{$sample2}{$marker} =~ /\.\/\./))
			{
				#print $variants{$sample1}{$marker}, ", $variants{$sample2}{$marker} ";
				#print "bad ";
				$misscount = $misscount + 1;
			}
			elsif ($variants{$sample1}{$marker} eq $variants{$sample2}{$marker})
			{
				#print "good ";
				$matchcount = $matchcount + 1;
				$totalcount = $totalcount + 1;
			}	
			else
			{
				#print "meh ";
				$totalcount = $totalcount + 1;
			}
			#print $variants{$sample1}{$marker}, ",$variants{$sample2}{$marker},$matchcount,$totalcount,$misscount ";
		}

		$percent = ($matchcount/$totalcount)*100;
		$subpercent = substr($percent,0,index($percent,'.') + 1 + 2);
		
		#begin printing the matrix!
		#for the full square matrix, uncomment the last command in the first if clause
		#for the half matrix, leave the last command in
		if ($sample1 eq $sample2)
		{
				if ($options{c} eq "off")
				{
					print OUTFILE "\t$matchcount", "/", $totalcount;
				}
				elsif ($options{c} eq "on")
				{
					print OUTFILE "\t$percent";
				}
				else
				{
					die "ERROR: invalid value for -c! Please run the script again with -c on or -c off.\n";
				}
				#last;
		}
		else
		{
			if ($options{c} eq "off")
			{
				print OUTFILE "\t$matchcount", "/", $totalcount;
			}
			elsif ($options{c} eq "on")
			{
				print OUTFILE "\t$percent";
			}
			else
			{
				die "ERROR: invalid value for -c! Please run the script again with -c on or -c off.\n";
			}
		}
	}
	#print "\t$count";
	print OUTFILE "\n";
}

#print Dumper($variants{"LR-26-55-B"});
