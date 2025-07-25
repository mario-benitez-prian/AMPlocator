import argparse
from amplocator.predictor import run_prediction

def main():

    parser = argparse.ArgumentParser(description="AMPlocator: two neural networks to predict antimicrobial peptide precursors and mature/active peptide location")
    parser.add_argument("--input_file", "-i", required = True, help="File in fasta format with protein data (amino acids)")
    parser.add_argument("--output_prefix", "-o", required = True, help="Output file name or path to save the results")
    parser.add_argument("--mode", "-m",
                        choices=["precursor", "locator", "full"], 
                        required=True,
                        help="Mode of prediction: 'precursors' (only precursor detection), 'locator' (only predict AMP mature peptide region in candidate precursors) or 'full' (both steps)")
    
    args = parser.parse_args()

    run_prediction(
        input_file=args.input_file,
        output_prefix=args.output_prefix,
        mode=args.mode
    )
