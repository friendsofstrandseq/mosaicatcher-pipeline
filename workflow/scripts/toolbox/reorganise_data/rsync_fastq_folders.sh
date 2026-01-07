# Script will rsync fastq folders from the mosaicatcher-pipeline to the scratch space for dedicated analysis
while IFS=$'\t' read -r date_folder subfolder; do
 date_folder=$(echo "$date_folder" | xargs)
 subfolder=$(echo "$subfolder" | xargs)
 src_path="/g/korbel/STOCKS_WF/mosaicatcher-pipeline/${date_folder}/${subfolder}"
 target_dir="/scratch/tweber/DATA/MC_DATA/EVA/EVA_MM10_250125/${subfolder}"
 mkdir -p "$target_dir"
 find "$src_path" -type d -name "fastq" | while read -r fastq_dir; do
   rsync -a "${fastq_dir}" "$target_dir"
 done
done <<EOF
2020-02-25-H2NCTAFX2	ondox5days20h20uMx02
2020-02-25-H2NCTAFX2	onhealthy20h20uMx02
2020-07-14-HCC35AFX2	PACO22PRnegx01
2020-07-14-HCC35AFX2	PACO22PRposx02
2021-02-17-HM7LYAFX2	OrgxCancerx02
2021-02-17-HM7LYAFX2	OrgxDoxocx02
2021-07-29-HWYJ5AFX2	10uMDXR42hx01
2021-07-29-HWYJ5AFX2	100uMDXR42hx01
2022-10-05-H2JJVAFX5	DXR100nMxS2x03
2022-09-26-H2LK5AFX5	DXR100nMxS2x01
2022-09-26-H2LK5AFX5	control30hS2x02
EOF