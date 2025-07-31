import pandas as pd
import numpy as np
import re
from amplocator.preprocess_data import preprocess_fasta_sequences
from amplocator.io_parser import read_fasta, write_fasta, write_precursor_predictions_table, write_locator_predictions_table, write_full_predictions_table

from pathlib import Path
import importlib.resources as pkg_resources
import amplocator  # Asegúrate que tu paquete se llama así

def get_model_path(filename):
    """Return the path to a model file inside the installed package."""

    with pkg_resources.as_file(pkg_resources.files(amplocator) / "cache/models" / filename) as path:
        return str(path)

def predict_precursors(headers, sequences, model_path):

    from tensorflow.keras.models import load_model

    print("\n----------RUNNING PRECURSOR MODEL----------\n")

    print("[INFO] Reading and preprocessing data for precursor prediction...")

    X = preprocess_fasta_sequences(sequences)

    model = load_model(model_path)

    print("[INFO] Predicting precursors...")
    preds = model.predict(X, verbose=1).flatten()

    results = pd.DataFrame({
        "ID": headers,
        "Precursor": sequences,
        "Precursor_score": preds
    })

    results = results[results["Precursor_score"] >= 0.5] # Filtering positive results

    num_positives = len(results["ID"])

    print(f"[INFO] Detected {num_positives} positive precursors.")

    results["Precursor_score"] = results["Precursor_score"].round(2) 

    return results




def predict_amp_regions(headers, sequences, model_path):

    from tensorflow.keras.models import load_model

    print("\n----------RUNNING LOCATOR MODEL----------\n")

    print("[INFO] Reading and preprocessing data for mature AMP localization...")
    X = preprocess_fasta_sequences(sequences)

    model = load_model(model_path)

    print("[INFO] Predicting mature AMP regions...")
    preds = model.predict(X, verbose=1)  # (N, L) probabilities

    results = []

    for i, pred in enumerate(preds):
        real_length = len(sequences[i])
        start_of_sequence = np.where(X[i].sum(axis=1) != 0)[0][0]
        pred_real = pred[start_of_sequence:start_of_sequence + real_length]
        amp_indices = np.where(pred_real >= 0.5)[0]

        '''print("\n--- Sequence:", headers[i])
        print("Original sequence:", sequences[i])
        print("Predicted full vector (length={}):".format(len(pred)), np.round(pred, 2).tolist())
        print("Trimmed prediction (no padding):", np.round(pred_real, 2).tolist())
        print("AMP indices:", amp_indices.tolist())'''

        if len(amp_indices) > 0:

            mature_amp = ''.join(sequences[i][idx] if idx in amp_indices else '-' for idx in range(real_length))
            avg_prob = np.mean(pred_real[amp_indices])  # Use pred_real instead of full pred
            '''print("Mature AMP sequence:", mature_amp)
            print("Mean AMP probability:", round(avg_prob, 2))'''
        else:
            mature_amp = '-' * real_length
            avg_prob = np.mean(pred_real)

        results.append({
            "ID": headers[i],
            "Precursor": sequences[i],
            "Mature_peptide": mature_amp,
            "Mature_score": round(avg_prob, 2)
        })

    results = pd.DataFrame(results)

    return results




def run_prediction(input_file, output_prefix, mode):

    import tensorflow as tf

    output_prefix = re.sub(r'\s*\(.*?\)', '', output_prefix)  # remove " (1)", etc.
    output_prefix = output_prefix.replace(" ", "_")           # replace spaces with "_"

    headers, sequences = read_fasta(input_file)

    if mode == "precursor":
        precursor_results = predict_precursors(headers, sequences, get_model_path("precursor_model.keras"))

        if precursor_results.empty:
            print("[WARNING] No precursors detected.")
            return
        else:
            write_fasta(precursor_results, output_prefix)
            write_precursor_predictions_table(precursor_results, output_prefix)

    elif mode == "full":
        precursor_results = predict_precursors(headers, sequences, get_model_path("precursor_model.keras"))

        if precursor_results.empty:
            print("[WARNING] No precursors detected.")
            return
        else:
            write_fasta(precursor_results, output_prefix)
            locator_results = predict_amp_regions(list(precursor_results["ID"]), list(precursor_results["Precursor"]), get_model_path("amp_locator_model.keras"))
            write_full_predictions_table(precursor_results, locator_results, output_prefix)

    elif mode == "locator":
        locator_results = predict_amp_regions(headers, sequences, get_model_path("amp_locator_model.keras"))

        if locator_results.empty:
            print("[WARNING] No mature AMP detected.")
            return
        else:
            write_locator_predictions_table(locator_results, output_prefix)

    print("[INFO] Prediction completed.")

