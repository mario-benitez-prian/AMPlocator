

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