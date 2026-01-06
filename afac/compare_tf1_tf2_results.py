import numpy as np
import tensorflow as tf

# Load predictions
tf1_prediction = np.load("tf1_prediction.npy", allow_pickle=True)
tf2_prediction = np.load("tf2_prediction.npy", allow_pickle=True)

# Debug: Inspect TF2 prediction type and content
print("Type of TF2 Prediction:", type(tf2_prediction))
print("TF2 Prediction Content:", tf2_prediction)

# Extract tensor from TF2 prediction (adjust key as necessary)
if isinstance(tf2_prediction, dict):
    tf2_prediction_array = tf2_prediction["dense_43"].numpy()
elif isinstance(tf2_prediction, tf.Tensor):
    tf2_prediction_array = tf2_prediction.numpy()
else:
    raise ValueError("Unexpected type for TF2 prediction.")

# Print the types and shapes of predictions
print(
    f"TF1 Prediction Type: {type(tf1_prediction)}, Shape: {getattr(tf1_prediction, 'shape', None)}"
)
print(
    f"TF2 Prediction Type: {type(tf2_prediction)}, Shape: {getattr(tf2_prediction, 'shape', None)}"
)

# Check the data in the predictions
print("TF1 Prediction Data:", tf1_prediction)
print("TF2 Prediction Data:", tf2_prediction)

# Compare predictions
# np.testing.assert_allclose(tf1_prediction, tf2_prediction, atol=1e-5)
# print("TF1 and TF2 predictions match within tolerance.")

# Convert predictions to NumPy arrays if needed
tf1_prediction = np.array(tf1_prediction)
tf2_prediction = np.array(tf2_prediction)

# Check for invalid values
if not np.isfinite(tf1_prediction).all():
    print("TF1 Prediction contains invalid values (NaN, inf, or None).")

if not np.isfinite(tf2_prediction).all():
    print("TF2 Prediction contains invalid values (NaN, inf, or None).")

# Compare predictions
np.testing.assert_allclose(tf1_prediction, tf2_prediction, atol=1e-5)
print("TF1 and TF2 predictions match within tolerance.")
