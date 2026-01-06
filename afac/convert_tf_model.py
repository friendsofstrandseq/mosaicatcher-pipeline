import numpy as np
import tensorflow as tf
from tensorflow.keras.models import load_model

# Path to the legacy H5 model
legacy_model_path = "/g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq/workflow/data/scNOVA/models_CNN/DNN_train80_chr22.h5"

# Load the legacy model without the optimizer state
model = load_model(legacy_model_path, compile=False)  # Prevent loading optimizer state


# Path to save the model in TensorFlow 2.x format
converted_model_path = "/g/korbel2/weber/workspace/StrandSeq_workspace/DEV/mosaicatcher-pipeline-friendsofstrandseq/workflow/data/scNOVA/models_CNN/DNN_train80_chr22_converted"

# Save the model without the optimizer
model.save(converted_model_path, save_format="tf")


print(f"Model input shape: {model.input_shape}")
print(f"Model successfully stripped of optimizer and saved to: {converted_model_path}")
# Generate and save input data
dummy_input = np.random.random((1, 150, 5)).astype(np.float32)
np.save("input_data.npy", dummy_input)
print("Input data saved.")

# Perform inference
tf1_prediction = model.predict(dummy_input)

# Save the prediction
np.save("tf1_prediction.npy", tf1_prediction)
print("TF1 prediction saved.")
