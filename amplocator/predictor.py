import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from amplocator.preprocess_data import preprocess_fasta_sequences
from amplocator.fasta_parser import read_fasta, write_fasta, write_precursor_predictions_table

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
    preds = model.predict(X, verbose=1)  # (N, L) probabilities

    binary_labels = (preds > 0.5).astype(int)

    mature_amps = []
    trimmed_labels = []
    trimmed_preds = []

    for seq, label, prob in zip(sequences, binary_labels, preds):
        length = len(seq)
        label = label[:length]
        prob = prob[:length]

        masked_seq = ''.join([aa if l == 1 else '-' for aa, l in zip(seq, label)])
        mature_amps.append(masked_seq)
        trimmed_labels.append(label)
        trimmed_preds.append(prob)

    return sequences, mature_amps, trimmed_labels, trimmed_preds

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
        full_seqs, mature_seqs, labels, scores = predict_amp_regions(sequences, max_length, locator_model_path)
        write_locator_predictions_table(headers, full_seqs, mature_seqs, labels, scores, output_prefix)
        write_fasta(headers, mature_seqs, output_prefix)

    print("[INFO] Prediction completed.")


    ########### Parte de lectura y escritura de archivos 

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    headers = [record.description for record in records]
    sequences = [str(record.seq).upper().strip() for record in records]
    return headers, sequences

def write_fasta(headers, sequences, output_prefix):
    records = [
        SeqRecord(Seq(seq), id=header, description="")
        for header, seq in zip(headers, sequences)
    ]
    fasta_file = f"{output_prefix}_precursor.fasta"
    SeqIO.write(records, f"{output_prefix}_precursor.fasta", "fasta")

def write_precursor_predictions_table(headers, seqs, preds, output_prefix):
    """
    Export results of precursor proteins
    """
    tsv_file = f"{output_prefix}_precursor.tsv"
    xlsx_file = f"{output_prefix}_precursor.xlsx"
    #print(headers)
    #print(seqs)
    #print(preds)
    #print(len(preds))
    #assert len(headers) == len(seqs) == len(preds), "Las longitudes no coinciden"

    df = pd.DataFrame({
        "ID": headers,
        "Precursor": seqs, 
        "Probability": preds
    }) 

    print(f"[INFO] Precursor results save in: {tsv_file} y {xlsx_file}")
    df.to_csv(tsv_file, sep="\t", index=False)
    df.to_excel(xlsx_file, index=False)

def write_locator_predictions_table(headers, full_seqs, mature_seqs, labels, scores, output_prefix):
    """
    Export results of mature AMP localization.

    """
    tsv_file = f"{output_prefix}_locator_predictions.tsv"
    xlsx_file = f"{output_prefix}_locator_predictions.xlsx"

    mean_scores = [np.mean(s) for s in scores] # Calculating mean probability 

    df = pd.DataFrame({
        "ID": headers,
        "Full_Seq": full_seqs,
        "Mature_AMP": mature_seqs,
        "Probability": mean_scores
    })

    df = df[df["Probability"] >= 0.5] # Filtering positive predictions only

    print(f"[INFO] Saving mature AMP prediction results to: {tsv_file} and {xlsx_file}")
    df.to_csv(tsv_file, sep="\t", index=False)
    df.to_excel(xlsx_file, index=False)

def write_full_predictions_table():
    """
    Export precursor and mature AMP results 

    """
    print("[INFO] Saving precursors and mature AMPs")