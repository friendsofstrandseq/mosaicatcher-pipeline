import pandas as pd

df = pd.read_csv(
    "{output_folder}/counts/{sample}/{sample}.info".format(
        output_folder="/g/korbel2/weber/workspace/mosaicatcher-update/.tests/output_CHR21_selected", sample="RPE-BM510"
    ),
    skiprows=13,
    sep="\t",
)
print(df)
cell_list = df.cell.tolist()
print(cell_list)

print(
    [
        sub_e
        for e in [
            expand(
                "{output_folder}/segmentation/{sample}/segmentation-per-cell/{cell}.txt",
                output_folder=config["output_location"],
                sample=sample,
                cell=cell_list,
            )
            for sample in samples
        ]
        for sub_e in e
    ]
)
