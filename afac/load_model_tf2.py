import tensorflow as tf
from keras.layers import TFSMLayer

# Path to the TensorFlow SavedModel directory
converted_model_path = "/g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq/workflow/data/scNOVA/models_CNN/DNN_train80_chr22_converted"

# Load the model using TFSMLayer
model = TFSMLayer(converted_model_path, call_endpoint="serving_default")

# Use the model for inference (example with dummy input)
import numpy as np

# Create dummy input matching the model's input shape
# Example for a Conv1D model: (batch_size, sequence_length, num_channels)
# Load input data
dummy_input = np.load("input_data.npy")
print("Input data loaded.")


# Perform inference
tf2_prediction = model(dummy_input)

# Save the prediction
np.save("tf2_prediction.npy", tf2_prediction)
print("TF2 prediction saved.")
