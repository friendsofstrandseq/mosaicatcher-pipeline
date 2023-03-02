#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from tensorflow import keras

input_args = sys.argv

Nucleosome_data_sample = pd.read_csv(input_args[1], delimiter="\t", header=None)
Nucleosome_data_new_all = Nucleosome_data_sample
TSS_matrix_new_all = pd.read_csv(input_args[2], delimiter="\t", header=None)

Nucleosome_data_new = Nucleosome_data_new_all
TSS_matrix_new = TSS_matrix_new_all


img_rows_new, img_cols_new = 150, 5
num_classes_new = 2
num_images_new = Nucleosome_data_new.shape[0]

x_as_array_new = Nucleosome_data_new.values[:, :]
x_shaped_array_new = x_as_array_new.reshape(num_images_new, img_rows_new, img_cols_new)
x_new = x_shaped_array_new
test_chr_new = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
]
test_chr_1st = ["chr1"]
test_chr_rest = [
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
]

##1) Load the model and perform prediction for each chromosomes
from tensorflow import keras
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Flatten, Conv2D, Conv1D, BatchNormalization, MaxPooling1D, GlobalMaxPooling1D, Dropout
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from pandas import DataFrame


filename = input_args[3]

i = input_args[4]

j = input_args[5]

# for i in test_chr_1st:
#     x_test = x_new[TSS_matrix_new[0] == i, :]
#     modelname = "utils/models_CNN/DNN_train{}_".format(j) + str(i) + ".h5"
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns=["prob1", "prob2"])
#     df_new = df

# for i in test_chr_rest:
x_test = x_new[TSS_matrix_new[0] == i, :]
modelname = "workflow/data/scNOVA/models_CNN/DNN_train{}_".format(j) + str(i) + ".h5"
Nucleosome_model_fixed = load_model(modelname)
y_pred = Nucleosome_model_fixed.predict_proba(x_test)
df = DataFrame(y_pred, columns=["prob1", "prob2"])

# df_new = df_new.append(df)
# for i in test_chr_1st:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train80_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df

# for i in test_chr_rest:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train80_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df_new.append(df)


df.to_csv(filename)


# for i in test_chr_1st:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train40_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df

# for i in test_chr_rest:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train40_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df_new.append(df)

# filename = input_args[4]
# df_new.to_csv(filename)


# for i in test_chr_1st:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train20_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df

# for i in test_chr_rest:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train20_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df_new.append(df)

# filename = input_args[5]
# df_new.to_csv(filename)


# for i in test_chr_1st:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train5_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df

# for i in test_chr_rest:
#     x_test = x_new[TSS_matrix_new[0]==i,:]
#     modelname = 'utils/models_CNN/DNN_train5_' + str(i) + '.h5'
#     Nucleosome_model_fixed = load_model(modelname)
#     y_pred = Nucleosome_model_fixed.predict_proba(x_test)
#     df = DataFrame(y_pred, columns= ['prob1', 'prob2'])
#     df_new = df_new.append(df)

# filename = input_args[6]
# df_new.to_csv(filename)
