REF="/g/korbel2/weber/MosaiCatcher_files/refgenomes_human_local/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
mkdir -p LITE/RPE-BM510/all/ LITE/RPE-BM510/selected/ TMP/
for file in *.bam; do
	echo "$file" &&
		cell=$(echo "${file}" | grep -P -o "PE[0-9]*") &&
		samtools view "$file" chr21 -b -o TMP/"$cell".bam &&
		bamToFastq -i TMP/"$cell".bam -fq TMP/"$cell".1.fq -fq2 TMP/"$cell".2.fq &&
		bwa mem -R "@RG\\tID:$cell\\tPL:Illumina\\tSM:RPE-BM510" "$REF" "${file%.*}".fq | samtools sort | samtools view -b -o LITE/RPE-BM510/all/"$file" &&
		cd LITE/RPE-BM510/selected/ &&
		ln -s ../all/* .
done
