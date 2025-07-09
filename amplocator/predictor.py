import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from amplocator.preprocess_data import preprocess_fasta_sequences
from amplocator.fasta_parser import read_fasta, write_fasta, export_predictions_table

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

def predict_amp_regions(sequences, max_length, model_path):
    print("[INFO] Preprocessing data for AMP localization...")

    X = preprocess_fasta_sequences(sequences, max_length)

    print("[INFO] Loading AMP locator model...")
    model = load_model(model_path)

    print("[INFO] Predicting AMP regions...")
    preds = model.predict(X, verbose=1)

    result = []
    for seq, prob in zip(sequences, preds):
        binary = (prob > 0.5).astype(int)
        if np.sum(binary) == 0:
            result.append("")
        else:
            start = np.argmax(binary)
            end = len(binary) - np.argmax(binary[::-1])
            result.append(seq[start:end])
    return result

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
        export_predictions_table(filtered_headers, filtered_seqs, filtered_preds, output_prefix)

    elif mode == "full":
        labels = predict_precursors(sequences, max_length, precursor_model_path)
        prec_headers = [h for h, l in zip(headers, labels) if l == 1]
        prec_seqs = [s for s, l in zip(sequences, labels) if l == 1]

        amp_seqs = predict_amp_regions(prec_seqs, max_length, locator_model_path)
        final_headers = [h for h, amp in zip(prec_headers, amp_seqs) if amp]
        final_seqs = [amp for amp in amp_seqs if amp]
        write_fasta(final_headers, final_seqs, output_file)

    elif mode == "locator":
        amp_seqs = predict_amp_regions(sequences, max_length, locator_model_path)
        final_headers = [h for h, amp in zip(headers, amp_seqs) if amp]
        final_seqs = [amp for amp in amp_seqs if amp]
        write_fasta(final_headers, final_seqs, output_file)

    print("[INFO] Prediction completed.")