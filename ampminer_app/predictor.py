import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from ampminer_app.preprocess_data import preprocess_fasta_sequences
from ampminer_app.fasta_parser import read_fasta, write_fasta

def predict_precursors(sequences, max_length, model_path):
    print("[INFO] Preprocessing data for precursor prediction...")

    X = preprocess_fasta_sequences(sequences, max_length)

    print(f"[INFO] Input shape: {X.shape}")

    print("[INFO] Loading precursor model...")
    model = load_model(model_path)

    print("[INFO] Predicting precursors...")
    preds = model.predict(X, verbose=0)
    labels = (preds > 0.5).flatten()
    print(f"[INFO] Detected {np.sum(labels)} positive precursors.")
    return labels

def predict_amp_regions(sequences, max_length, model_path):
    print("[INFO] Preprocessing data for AMP localization...")

    X = preprocess_fasta_sequences(sequences, max_length)

    print("[INFO] Loading AMP locator model...")
    model = load_model(model_path)

    print("[INFO] Predicting AMP regions...")
    preds = model.predict(X, verbose=0)

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

def run_prediction(fasta_file, output_file, mode):
    print("[INFO] Reading input FASTA...")
    headers, sequences = read_fasta(fasta_file)

    max_length = 300
    precursor_model_path = "models/precursor_model.keras"
    locator_model_path = "models/amp_locator_model.keras"

    if mode == "precursor":
        labels = predict_precursors(sequences, max_length, precursor_model_path)
        filtered_headers = [h for h, l in zip(headers, labels) if l == 1]
        filtered_seqs = [s for s, l in zip(sequences, labels) if l == 1]
        write_fasta(filtered_headers, filtered_seqs, output_file)

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
