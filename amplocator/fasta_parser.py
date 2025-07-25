

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq




def read_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    headers = [record.description for record in records]
    sequences = [str(record.seq).upper().strip() for record in records]
    return headers, sequences




def write_fasta(results, output_prefix):
    records = [
        SeqRecord(Seq(seq), id=header, description="")
        for header, seq in zip(results["ID"], results["Precursor"])
    ]
    fasta_file = f"{output_prefix}_precursor.fasta"
    SeqIO.write(records, f"{output_prefix}_precursor.fasta", "fasta")




def write_precursor_predictions_table(results, output_prefix):
    """
    Export results of precursor proteins

    """
    if results.empty:
        print("[WARNING] No mature AMP predictions to export.")
        return

    tsv_file = f"{output_prefix}_precursor.tsv"
    xlsx_file = f"{output_prefix}_precursor.xlsx"

    print(f"[INFO] Precursor results save in: {tsv_file} y {xlsx_file}")
    results.to_csv(tsv_file, sep="\t", index=False)
    results.to_excel(xlsx_file, index=False)




def write_locator_predictions_table(results, output_prefix):
    """
    Export results of mature AMP localization from a unified results list.

    """
    if results.empty:
        print("[WARNING] No mature AMP predictions to export.")
        return

    tsv_file = f"{output_prefix}_locator_predictions.tsv"
    xlsx_file = f"{output_prefix}_locator_predictions.xlsx"

    print(f"[INFO] Saving mature AMP prediction results to: {tsv_file} and {xlsx_file}")
    results.to_csv(tsv_file, sep="\t", index=False)
    results.to_excel(xlsx_file, index=False)




def write_full_predictions_table(precursor_results, locator_results, output_prefix):
    """
    Export precursor and mature AMP results 

    """
    
    if precursor_results.empty:
        print("[WARNING] No mature AMP predictions to export.")
        return

    tsv_file = f"{output_prefix}_full_predictions.tsv"
    xlsx_file = f"{output_prefix}_full_predictions.xlsx"

    # Seleccionamos solo las columnas que queremos del locator_results
    locator_subset = locator_results[["Mature_peptide", "Mature_score"]].reset_index(drop=True)

    # Nos aseguramos de que precursor_results también tenga índice limpio
    precursor_results = precursor_results.reset_index(drop=True)

    # Concatenamos columna a columna (axis=1)
    full_results = pd.concat([precursor_results, locator_subset], axis=1)

    '''full_results = pd.DataFrame({
        "ID": precursor_results["ID"],
        "Precursor": precursor_results["Precursor"],
        "Precursor_score": precursor_results["Precursor_score"],
        "Mature_peptide": locator_results["Mature_peptide"],
        "Mature_score": locator_results["Mature_score"]
    })'''
    print(full_results.head(10))
    print(full_results.shape)

    print(f"[INFO] Saving precursors and mature AMPs results to: {tsv_file} and {xlsx_file}")
    full_results.to_csv(tsv_file, sep="\t", index=False)
    full_results.to_excel(xlsx_file, index=False)

    print(full_results.head(10))
    print(full_results.shape)