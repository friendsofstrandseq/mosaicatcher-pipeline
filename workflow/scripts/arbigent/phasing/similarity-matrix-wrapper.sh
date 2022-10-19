#!/bin/bash
while getopts c:f:o: flag
	do
	    case "${flag}" in
	        c) calc_in_percent=${OPTARG};;
	        f) vcffile=${OPTARG};;
	        o) outfile=${OPTARG};;
	    esac
	done
	echo "calc_in_percent: $calc_in_percent";
	echo "vcffile: $vcffile";
	echo "outfile: $outfile";

nsnps=$(bcftools view -H "$vcffile" | head | wc -l)

if ((nsnps > 0)); then
    echo "ok im perling now"
	/usr/bin/perl scripts_phasing/similarity-matrix.pl \
		-c "$calc_in_percent" \
		-f "$vcffile" \
		-o $outfile"_noname"
	sed ':a;N;$!ba;s/\n/	/g' <(cat <(echo $vcffile) $outfile"_noname") > $outfile

else
    echo "no this is so small"
	echo -e "$vcffile	sample_name	none	none	none	100	0	none	0	100" > $outfile
fi

echo "$nsnps"
