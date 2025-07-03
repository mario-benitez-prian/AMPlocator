

import argparse
from ampminer_app.predictor import run_prediction

def main():
    parser = argparse.ArgumentParser(description="Predicción de precursores AMP en proteoma completo.")
    parser.add_argument("fasta_file", help="Archivo FASTA de entrada con el proteoma completo.")
    parser.add_argument("model_file", help="Ruta al modelo entrenado (formato .h5).")
    parser.add_argument("output_file", help="Archivo FASTA de salida con las secuencias positivas.")
    
    args = parser.parse_args()

    run_prediction(
        fasta_file=args.fasta_file,
        model_file=args.model_file,
        output_file=args.output_file
    )
    
print("HOLA SOY CLI")
