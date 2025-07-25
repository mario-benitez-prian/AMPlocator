import argparse
from amplocator.predictor import run_prediction

def main():

    parser = argparse.ArgumentParser(description="AMPlocator: two neural networks to predict antimicrobial peptide precursors and mature/active peptide location")
    parser.add_argument("fasta_file", help="File in fasta format with protein data")
    parser.add_argument("output_prefix", help="Output fasta file with predicted sequences")
    parser.add_argument("--mode", 
                        choices=["precursor", "full", "locator"], 
                        required=True,
                        help="Mode of prediction: 'precursors' (only precursor detection), 'full' (both steps), or 'locator' (AMP localization only)")
    
    args = parser.parse_args()

    run_prediction(
        fasta_file=args.fasta_file,
        output_prefix=args.output_prefix,
        mode=args.mode
    )
