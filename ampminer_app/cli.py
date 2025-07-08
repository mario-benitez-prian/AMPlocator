import argparse
from ampminer_app.predictor import run_prediction

def main():
    parser = argparse.ArgumentParser(description="Predicción de precursores AMP en proteoma completo.")
    parser.add_argument("fasta_file", help="Archivo FASTA de entrada con el proteoma completo.")
    parser.add_argument("output_file", help="Archivo FASTA de salida con las secuencias positivas.")
    parser.add_argument("--mode", 
                        choices=["precursor", "full", "locator"], 
                        required=True,
                        help="Mode of prediction: 'precursors' (only precursor detection), 'full' (both steps), or 'locator' (AMP localization only)")
    
    args = parser.parse_args()

    run_prediction(
        fasta_file=args.fasta_file,
        output_file=args.output_file,
        mode=args.mode
    )
    
print("HOLA SOY CLI")
