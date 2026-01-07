import pandas as pd

l = list()
for file in list(sorted(snakemake.input.bam)):
    tmp_idxstats = pd.read_csv(file, sep="\t", names=["chrom", "seq_len", "reads_mapped", "reads_unmapped"], header=False)
    l.append({"Sample": snakemake.wildcards.sample, "Cell": snakemake.wildcards.cell, "Total_reads": int(tmp_idxstats.reads_mapped.sum())})
final_df = pd.concat(df).sort_values(by="cell")
final_df.to_csv(snakemake.output[0])
