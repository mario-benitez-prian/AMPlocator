

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    headers = [record.description for record in records]
    sequences = [str(record.seq) for record in records]
    return headers, sequences

def write_fasta(headers, sequences, output_prefix):
    records = [
        SeqRecord(Seq(seq), id=header, description="")
        for header, seq in zip(headers, sequences)
    ]
    fasta_file = f"{output_prefix}_precursor.fasta"
    SeqIO.write(records, f"{output_prefix}_precursor.fasta", "fasta")

def export_predictions_table(filtered_headers, filtered_seqs, filtered_preds, output_prefix):
    
    tsv_file = f"{output_prefix}_precursor.tsv"
    xlsx_file = f"{output_prefix}_precursor.xlsx"
    #print(filtered_headers)
    #print(filtered_seqs)
    #print(filtered_preds)
    #print(len(filtered_preds))
    #assert len(filtered_headers) == len(filtered_seqs) == len(filtered_preds), "Las longitudes no coinciden"

    df = pd.DataFrame({
        "ID": filtered_headers,
        "Precursor": filtered_seqs, 
        "Probability": filtered_preds
    })

    print(f"[INFO] Precursor results save in: {tsv_file} y {xlsx_file}")
    df.to_csv(tsv_file, sep="\t", index=False)
    df.to_excel(xlsx_file, index=False)