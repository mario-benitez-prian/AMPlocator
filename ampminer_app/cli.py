import argparse
from ampminer_app.predictor import run_prediction

def main():
    parser = argparse.ArgumentParser(description="Predicción de precursores AMP en proteoma completo.")
    parser.add_argument("fasta_file", help="Archivo FASTA de entrada con el proteoma completo.")
    parser.add_argument("model_file", help="Ruta al modelo entrenado (formato .h5).")
    parser.add_argument("output_file", help="Archivo FASTA de salida con las secuencias positivas.")
    parser.add_argument("--max_length", type=int, default=300, help="Longitud máxima para las secuencias (por defecto: 300).")
    
    args = parser.parse_args()

    run_prediction(
        fasta_file=args.fasta_file,
        model_file=args.model_file,
        output_file=args.output_file,
        max_length=args.max_length
    )
