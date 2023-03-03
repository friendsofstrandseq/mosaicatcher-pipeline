import pandas as pd
from tensorflow.keras.models import load_model

Nucleosome_data_new = pd.read_csv(
    snakemake.input.features,
    delimiter="\t",
    header=None,
)
TSS_matrix_new = pd.read_csv(
    snakemake.input.TSS_annot,
    delimiter="\t",
    header=None,
)

i = str(snakemake.wildcards.chrom)

j = str(snakemake.wildcards.i)

TSS_matrix_new_index = TSS_matrix_new.loc[TSS_matrix_new[0] == i].index.tolist()
x_test = Nucleosome_data_new.loc[TSS_matrix_new_index].values.reshape(len(TSS_matrix_new_index), 150, 5)

Nucleosome_model_fixed = load_model("workflow/data/scNOVA/models_CNN/DNN_train{}_".format(j) + str(i) + ".h5")
y_pred = Nucleosome_model_fixed.predict_proba(x_test)
df = pd.DataFrame(y_pred, columns=["prob1", "prob2"])

df.to_csv(snakemake.output.train)
