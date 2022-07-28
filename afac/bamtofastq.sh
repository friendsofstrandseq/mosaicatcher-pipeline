dir=$1
for file in "$dir"/*.sort.mdup.bam; do
  echo "$file"
  # samtools view -o "${file%.*}".chr21.bam "${file}" chr21
  # samtools sort -n -o "${file%.*}".chr21.sort.bam "${file%.*}".chr21.bam
  # bamToFastq -i "${file%.*}".chr21.sort.bam -fq "${file%.*}".1.fastq -fq2 "${file%.*}".2.fastq
  # bgzip "${file%.*}".1.fastq && bgzip "${file%.*}".2.fastq
  # rm "${file%.*}".chr21.bam "${file%.*}".chr21.sort.bam
  echo "${file%%.*}"
  mv "${file%.*}".1.fastq.gz "${file%%.*}".1.fastq.gz
  mv "${file%.*}".2.fastq.gz "${file%%.*}".2.fastq.gz

done
