

import tensorflow as tf
from ampminer_app.preprocess_data import preprocess_fasta_sequences
from ampminer_app.fasta_parser import load_fasta_as_dataframe, write_fasta

#tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

def run_prediction(fasta_file, model_file, output_file, max_length):
    print("[INFO] Loading model...")
    model = tf.keras.models.load_model(model_file)

    print("[INFO] Loading sequences...")
    df = load_fasta_as_dataframe(fasta_file)

    print("[INFO] Preprocessing sequences...")
    X = preprocess_data(df, max_length)

    print("[INFO] Predicting...")
    predictions = model.predict(X)
    predicted_labels = (predictions > 0.5).astype(int).flatten()

    print("[INFO] Filtering positive sequences...")
    df["Predicted"] = predicted_labels
    df_positive = df[df["Predicted"] == 1]

    print(f"[INFO] Saving {len(df_positive)} results in: {output_file}")
    write_fasta(df_positive, output_file)

    print("[INFO] Work done!.")
