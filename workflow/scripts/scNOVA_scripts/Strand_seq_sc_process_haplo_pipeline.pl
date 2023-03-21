@files = <../bam/*.bam>;

foreach $file (@files) {
`samtools view -H $file > ../bam_haplo/header_test.sam`;
`samtools view -f 99 $file | cat ../bam_haplo/header_test.sam - | samtools view -Sb - > ../bam_haplo/$file.C1.bam`;
`samtools view -f 147 $file | cat ../bam_haplo/header_test.sam - | samtools view -Sb - > ../bam_haplo/$file.C2.bam`;
`samtools view -f 83 $file | cat ../bam_haplo/header_test.sam - | samtools view -Sb - > ../bam_haplo/$file.W1.bam`;
`samtools view -f 163 $file | cat ../bam_haplo/header_test.sam - | samtools view -Sb - > ../bam_haplo/$file.W2.bam`;
`samtools merge ../bam_haplo/$file.C.bam ../bam_haplo/$file.C1.bam ../bam_haplo/$file.C2.bam`;
`samtools merge ../bam_haplo/$file.W.bam ../bam_haplo/$file.W1.bam ../bam_haplo/$file.W2.bam`;
`samtools index ../bam_haplo/$file.C.bam`;
`samtools index ../bam_haplo/$file.W.bam`;
`rm ../bam_haplo/$file.C1.bam`;
`rm ../bam_haplo/$file.C2.bam`;
`rm ../bam_haplo/$file.W1.bam`;
`rm ../bam_haplo/$file.W2.bam`;

}
