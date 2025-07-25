import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from amplocator.preprocess_data import preprocess_fasta_sequences
from amplocator.fasta_parser import read_fasta, write_fasta, write_precursor_predictions_table, write_locator_predictions_table, write_full_predictions_table

import sys 
import os 
import logging
import tensorflow as tf

def predict_precursors(headers, sequences, max_length, model_path):

    print("[INFO] Preprocessing data for precursor prediction...")

    X = preprocess_fasta_sequences(sequences, max_length)

    model = load_model(model_path)

    print("[INFO] Predicting precursors...")
    preds = model.predict(X, verbose=1).flatten()
    
    labels = (preds >= 0.5).flatten()
    print(f"[INFO] Detected {np.sum(labels)} positive precursors.")

    results = pd.DataFrame({
        "ID": headers,
        "Precursor": sequences,
        "Precursor_score": preds.round(2)
    })

    return results

def predict_amp_regions(headers, sequences, max_length, model_path):

    print("[INFO] Preprocessing data for AMP localization...")
    X = preprocess_fasta_sequences(sequences, max_length)

    model = load_model(model_path)

    print("[INFO] Predicting AMP regions...")
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


def run_prediction(fasta_file, output_prefix, mode):

    print("[INFO] Reading input FASTA...")
    headers, sequences = read_fasta(fasta_file)

    max_length = 300
    precursor_model_path = "models/precursor_model.keras"
    locator_model_path = "models/amp_locator_model.keras"

    if mode == "precursor":
        results = predict_precursors(headers, sequences, max_length, precursor_model_path)
        results = results[results["Precursor_score"] >= 0.5] # Filter positive sequences
        write_fasta(results, output_prefix)
        write_precursor_predictions_table(results, output_prefix)

    elif mode == "full":
        precursor_results = predict_precursors(headers, sequences, max_length, precursor_model_path)
        precursor_results = precursor_results[precursor_results["Precursor_score"] >= 0.5] # Filter positive sequences

        write_fasta(precursor_results, output_prefix)

        locator_results = predict_amp_regions(list(precursor_results["ID"]), list(precursor_results["Precursor"]), max_length, locator_model_path)

        write_full_predictions_table(precursor_results, locator_results, output_prefix)

    elif mode == "locator":
        results = predict_amp_regions(headers, sequences, max_length, locator_model_path)
        write_locator_predictions_table(results, output_prefix)

    print("[INFO] Prediction completed.")

