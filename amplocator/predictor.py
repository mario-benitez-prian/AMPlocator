import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from amplocator.preprocess_data import preprocess_fasta_sequences
from amplocator.fasta_parser import read_fasta, write_fasta, write_precursor_predictions_table, write_locator_predictions_table

import sys 
import os 
import logging
import tensorflow as tf

def predict_precursors(sequences, max_length, model_path):

    print("[INFO] Preprocessing data for precursor prediction...")

    X = preprocess_fasta_sequences(sequences, max_length)

    print("[INFO] Loading precursor model...")
    model = load_model(model_path)

    print("[INFO] Predicting precursors...")
    preds = model.predict(X, verbose=1)
    labels = (preds > 0.5).flatten()

    print(f"[INFO] Detected {np.sum(labels)} positive precursors.")

    # Flatten preds to be unidimensional array
    preds = preds.flatten()

    return labels, preds

def predict_amp_regions(headers, sequences, max_length, model_path):
    print("[INFO] Preprocessing data for AMP localization...")
    X = preprocess_fasta_sequences(sequences, max_length)

    print("[INFO] Loading AMP locator model...")
    model = load_model(model_path)

    print("[INFO] Predicting AMP regions...")
    preds = model.predict(X, verbose=1)  # (N, L) probabilities

    results = []

    for i, pred in enumerate(preds):
        real_length = len(sequences[i])
        start_of_sequence = np.where(X[i].sum(axis=1) != 0)[0][0]
        pred_real = pred[start_of_sequence:start_of_sequence + real_length]
        amp_indices = np.where(pred_real > 0.5)[0]

        print("\n--- Sequence:", headers[i])
        print("Original sequence:", sequences[i])
        print("Predicted full vector (length={}):".format(len(pred)), np.round(pred, 2).tolist())
        print("Trimmed prediction (no padding):", np.round(pred_real, 2).tolist())
        print("AMP indices:", amp_indices.tolist())

        if len(amp_indices) > 0:
            mature_amp = ''.join(sequences[i][idx] if idx in amp_indices else '-' for idx in range(real_length))
            avg_prob = np.mean(pred_real[amp_indices])  # Use pred_real instead of full pred
            print("Mature AMP sequence:", mature_amp)
            print("Mean AMP probability:", round(avg_prob, 2))
            results.append((headers[i], sequences[i], mature_amp, round(avg_prob, 2)))
        else:
            print("No AMP detected.")

    print("\n[INFO] Final results:")
    #for r in results:
        #print(r)

    return results


def run_prediction(fasta_file, output_prefix, mode):

    print("[INFO] Reading input FASTA...")
    headers, sequences = read_fasta(fasta_file)

    max_length = 300
    precursor_model_path = "models/precursor_model.keras"
    locator_model_path = "models/amp_locator_model.keras"

    if mode == "precursor":
        labels, preds = predict_precursors(sequences, max_length, precursor_model_path)
        # Filter positive sequences
        filtered_headers = [h for h, l in zip(headers, labels) if l == 1]
        filtered_seqs = [s for s, l in zip(sequences, labels) if l == 1]
        filtered_preds = [s for s, l in zip(preds, labels) if l == 1]
        write_fasta(filtered_headers, filtered_seqs, output_prefix)
        write_precursor_predictions_table(filtered_headers, filtered_seqs, filtered_preds, output_prefix)

    elif mode == "full":
        labels = predict_precursors(sequences, max_length, precursor_model_path)
        prec_headers = [h for h, l in zip(headers, labels) if l == 1]
        prec_seqs = [s for s, l in zip(sequences, labels) if l == 1]

        amp_seqs = predict_amp_regions(prec_seqs, max_length, locator_model_path)
        final_headers = [h for h, amp in zip(prec_headers, amp_seqs) if amp]
        final_seqs = [amp for amp in amp_seqs if amp]
        write_fasta(final_headers, final_seqs, output_file)

    elif mode == "locator":
        results = predict_amp_regions(headers, sequences, max_length, locator_model_path)
        write_locator_predictions_table(results, output_prefix)

    print("[INFO] Prediction completed.")

