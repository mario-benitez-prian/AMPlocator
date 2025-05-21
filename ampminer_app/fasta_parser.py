import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def load_fasta_as_dataframe(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    df = pd.DataFrame({
        "ID": [record.description for record in records],
        "Full_Seq": [str(record.seq) for record in records]
    })
    return df

def write_fasta(df, output_file):
    records = [
        SeqRecord(Seq(row["Full_Seq"]), id=row["ID"], description="")
        for _, row in df.iterrows()
    ]
    SeqIO.write(records, output_file, "fasta")
