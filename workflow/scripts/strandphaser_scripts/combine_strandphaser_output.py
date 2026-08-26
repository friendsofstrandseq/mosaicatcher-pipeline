import pandas as pd

# Column schema of a per-chromosome phased_haps.txt, used only when there is nothing
# to combine (see the empty-input case below).
PHASED_HAPS_COLUMNS = [
    "sample",
    "cell",
    "chrom",
    "start",
    "end",
    "class",
    "hap1.cis.simil",
    "hap1.trans.simil",
    "hap2.cis.simil",
    "hap2.trans.simil",
]

l_df = list()
print(snakemake.input.files)
for j, file in enumerate(snakemake.input.files):
    print(j, file)
    tmp_df = pd.read_csv(file, sep="\t")
    print(tmp_df)
    l_df.append(tmp_df)

if l_df:
    df = pd.concat(l_df).drop_duplicates()
elif snakemake.config.get("allow_empty_phasing", False):
    # allow_empty_phasing=True: this sample has NO phaseable chromosomes (every
    # chromosome was filtered out by the min_snps_for_phasing gate — typically a very
    # low cell count -> too few SNVs). Emit an empty, header-only table so the sample
    # proceeds with empty phasing instead of crashing on pd.concat([]). Downstream
    # convert_strandphaser_output already tolerates an empty table, so the sample can
    # still reach arbigent (which genotypes predefined segments and does not require
    # phasing). Samples WITH phaseable chromosomes are unaffected (l_df is non-empty).
    df = pd.DataFrame(columns=PHASED_HAPS_COLUMNS)
else:
    # Default behaviour: fail with a clear, actionable message instead of the cryptic
    # "No objects to concatenate". Set allow_empty_phasing=True to let such samples
    # through with empty phasing.
    raise ValueError(
        "No phaseable chromosomes for this sample: every chromosome has fewer than "
        "min_snps_for_phasing SNVs (typically a very low cell count). Set "
        "allow_empty_phasing=True to proceed with empty phasing and still run arbigent."
    )

df.to_csv(snakemake.output[0], sep="\t", index=False)
