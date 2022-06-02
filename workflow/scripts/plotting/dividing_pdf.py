import pandas as pd

config_df = pd.read_csv(snakemake.params.config_df, sep="\t")

# tmp_dict = (
#     config_df.loc[config_df["all/selected"] == "selected", ["Sample", "Cell"]]
#     .groupby("Sample")["Cell"]
#     .apply(lambda r: sorted(list(r)))
#     .to_dict()
# )

tmp_dict = (
    config_df.loc[config_df["Selected"] == True, ["Sample", "Cell"]]
    .groupby("Sample")["Cell"]
    .apply(lambda r: sorted(list(r)))
    .to_dict()
)
tmp_dict = {s: {i + 1: c for i, c in enumerate(cell_list)} for s, cell_list in tmp_dict.items()}
for s in tmp_dict.keys():
    tmp_dict[s][0] = "SummaryPage"

import os
from PyPDF2 import PdfFileWriter, PdfFileReader

inputpdf = PdfFileReader(snakemake.input[0], "rb")

cell_name = tmp_dict[snakemake.wildcards.sample][int(snakemake.wildcards.i)]

output = PdfFileWriter()
output.addPage(inputpdf.getPage(int(snakemake.wildcards.i)))

tmp_output_path = os.path.dirname(snakemake.input[0]) + "/{}.{}.pdf".format(cell_name, snakemake.wildcards.i)

with open(tmp_output_path, "wb") as outputStream:
    output.write(outputStream)
