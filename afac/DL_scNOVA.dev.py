import pandas as pd
import tensorflow as tf


Nucleosome_data_new = pd.read_csv(
    "/g/korbel/Costea/Computational/StrandSeq/P6_Ini/PDX1609p6Ini/scNOVA_result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_clone1.txt",
    delimiter="\t",
    header=None,
)
TSS_matrix_new = pd.read_csv(
    "/g/korbel/Costea/Computational/StrandSeq/P6_Ini/PDX1609p6Ini/scNOVA_result/Features_reshape_all_TSS_matrix_woM_all_RT_clone1.txt",
    delimiter="\t",
    header=None,
)

i = str("chr22")

j = str("80")

TSS_matrix_new_index = TSS_matrix_new.loc[TSS_matrix_new[0] == i].index.tolist()
x_test = Nucleosome_data_new.loc[TSS_matrix_new_index].values.reshape(
    len(TSS_matrix_new_index), 150, 5
)

Nucleosome_model_fixed = tf.keras.models.load_model(
    "/g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq/workflow/data/scNOVA/models_CNN/DNN_train{}_".format(
        j
    )
    + str(i)
    + ".h5"
)
y_pred = Nucleosome_model_fixed.predict(x_test)
df = pd.DataFrame(y_pred, columns=["prob1", "prob2"])

df.to_csv("test.scNOVADL.infer.txt")
