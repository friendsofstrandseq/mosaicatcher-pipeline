# Need to be in the query directory where the initial bam are located
## Extract chr21 chromosome from bam file
## Reheader to keep only chr21 tag (used by mosaic count)
## Index new file


for file in *.bam; do;
    echo "$file" 
    samtools view "$file" chr21 -b > /g/korbel2/weber/MosaiCatcher_files/bam_KG_chr21/RPE1-WT/all/"$file"
    /home/tweber/.conda/envs/strandseqnation/bin/samtools reheader -c 'grep -P "^@HD|^@RG|^\@PG|^\@SQ\tSN:chr21"' /g/korbel2/weber/MosaiCatcher_files/bam_KG_chr21/RPE-BM510/selected/"$file" > /g/korbel2/weber/MosaiCatcher_files/bam_KG_chr21_RH/RPE-BM510/selected/"$file" 
    samtools index /g/korbel2/weber/MosaiCatcher_files/bam_KG_chr21_RH/RPE-BM510/selected/"$file"
done;