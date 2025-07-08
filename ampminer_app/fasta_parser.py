import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    headers = [record.description for record in records]
    sequences = [str(record.seq) for record in records]
    return headers, sequences

def write_fasta(headers, sequences, output_file):
    records = [
        SeqRecord(Seq(seq), id=header, description="")
        for header, seq in zip(headers, sequences)
    ]
    SeqIO.write(records, output_file, "fasta")

print(read_fasta("data/Arabidopsis_thaliana_rRNA_aminoacidos.fasta"))