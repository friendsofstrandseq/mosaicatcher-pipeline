#!/bin/bash

module load tabix/0.2.6-GCCcore-11.2.0


# hap_file = $1
vcf_file="$1"
echo "$vcf_file"
sample=$(echo "$vcf_file" | sed "s/strandphaser.*//g" | rev | cut -d'/' -f 2 | rev)

# touch "$vcf_file"

echo '##fileformat=VCFv4.2
##fileDate=2023-01-31
##source=StrandPhase_algorithm
##reference=BSgenome.Hsapiens.UCSC.hg38
##phasing=Strand-seq
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=Q1,Number=1,Type=Float,Description="Quality measure of allele 1 (1-entropy)*coverage">
##FORMAT=<ID=Q2,Number=1,Type=Float,Description="Quality measure of allele 2 (1-entropy)*coverage">
##FORMAT=<ID=P1,Number=1,Type=Float,Description="Probability value of allele 1">
##FORMAT=<ID=P2,Number=1,Type=Float,Description="Probability value of allele 2">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'$sample'
' | bgzip >"$vcf_file"
chmod 777 "$vcf_file"
