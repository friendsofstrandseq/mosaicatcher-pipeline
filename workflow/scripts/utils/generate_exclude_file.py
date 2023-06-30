import pandas as pd
import pysam

# import parmap
import multiprocessing as mp
import logging

# Set up logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")
# logging.basicConfig(filename="/Users/frank/mosaicatcher-pipeline/debug.log", level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")


m = mp.Manager()
l_df = m.list()
# l_df = list()


def filter_chrom(bam, l):
    # READ BAM FILE HEADER OF FIRST BAM IN THE PANDAS DF
    h = pysam.view("-H", bam)
    # h = pysam.view("-H", os.listdir(snakemake.input.bam + "selected")[0])
    h = [e.split("\t") for e in h.split("\n") if "@SQ" in e]

    l.extend(h)


parmap.starmap(filter_chrom, list(zip(list(snakemake.input.bam))), l_df, pm_pbar=False, pm_processes=10)
# for file in list(snakemake.input.bam):
#     filter_chrom(file, l_df)


# CONVERT TO PANDAS DF
df_h = pd.DataFrame((list(l_df)), columns=["TAG", "Contig", "LN"])

logging.info(f"Processed raw DataFrame: \n {df_h.to_string()}")


# PROCESS CONTIGS
output_h = pd.DataFrame(df_h["Contig"].str.replace("SN:", ""))


logging.info(f'List of chromosomes provided in the configuration: \n {snakemake.params["chroms"]}')


output_h = output_h.loc[~output_h["Contig"].isin(snakemake.params["chroms"])]


# Log the content of the output_h DataFrame
logging.info(f"Processed list of chromosomes to be removed: \n {output_h.to_string()}")


# EXPORT
output_h = output_h["Contig"].drop_duplicates()

logging.info(f"Processed list of chromosomes to be removed without duplicates: \n {output_h.to_string()}")


output_h.to_csv(snakemake.output[0], index=False, sep="\t", header=False)
