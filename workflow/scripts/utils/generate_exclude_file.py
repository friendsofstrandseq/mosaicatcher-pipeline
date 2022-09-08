import pandas as pd
import pysam
import os, sys
import parmap
import multiprocessing as mp

# df_config_files = pd.read_csv(, sep="\t")
m = mp.Manager()
l_df = m.list()

# bam_list = [
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/HG00268x01PE20401.sort.mdup.bam",
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/HG00512_I_015.sort.mdup.bam",
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/HG00514_IV_008.sort.mdup.bam",
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/HG01352x02PE20402.sort.mdup.bam",
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/NA12878_cell113.sort.mdup.bam",
#     "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all/NA19239_III_066.sort.mdup.bam",
# ]


def filter_chrom(bam, l):
    # READ BAM FILE HEADER OF FIRST BAM IN THE PANDAS DF
    h = pysam.view("-H", bam)
    # h = pysam.view("-H", os.listdir(snakemake.input.bam + "selected")[0])
    h = [e.split("\t") for e in h.split("\n") if "@SQ" in e]

    l.extend(h)


# parmap.starmap(filter_chrom, list(zip(list(bam_list))), l_df, pm_pbar=True, pm_processes=10)
parmap.starmap(filter_chrom, list(zip(list(snakemake.input.bam))), l_df, pm_pbar=False, pm_processes=10)

# CONVERT TO PANDAS DF
df_h = pd.DataFrame((list(l_df)), columns=["TAG", "Contig", "LN"])

# PROCESS CONTIGS
output_h = pd.DataFrame(df_h["Contig"].str.replace("SN:", ""))
# print(snakemake.params["chroms"], type(snakemake.params["chroms"]))
# chroms = ["chr{}".format(str(c)) for c in list(range(1, 23))] + ["chrX"]

output_h = output_h.loc[~output_h["Contig"].isin(snakemake.params["chroms"])]

# output_h = output_h.loc[~output_h["Contig"].isin(chroms)]

# print(output_h)

# EXPORT
output_h["Contig"].unique().to_csv(snakemake.output[0], index=False, sep="\t", header=False)
