# Make functions Bash command-line-friendly by parsing arguments
# - i.e., specify command-line arguments and define how they should be interpreted into python script
import argparse
import sys
sys.path.append('C:\\Users\\joanaq\\Documents\\learning-bioinformatics\\projects\\ros-orf')
from random_seq import random_dnaseq_generator, seqrecord_fasta_writer, seqrecord_printer

def main():
    parser = argparse.ArgumentParser(description="Generate random DNA sequences and either print them (default) or save to FASTA file (if file name specified).")
    parser.add_argument('-n', '--n_seq', type=int, required=True, help="Number of sequences to generate.")
    parser.add_argument('-l', '--length', type=int, required=True, help="Length of each sequence in bp.")
    parser.add_argument('-g', '--gc_pct', type=int, required=True, help="GC percentage (0-100).")
    parser.add_argument('-i', '--batch_id', type=str, required=True, help="General identifier for sequences.")
    parser.add_argument('-f', '--fasta', type=str, required=False, help="Output FASTA filename (omit to just print sequences).")

    args = parser.parse_args()

    # Generate random DNA sequences
    records = random_dnaseq_generator(args.n_seq, args.length, args.gc_pct, args.batch_id)

    # Print sequences or save them in FASTA file
    if args.fasta:
        seqrecord_fasta_writer(records, args.fasta)
        print(f"Generated {args.n_seq} random DNA sequences and saved to {args.fasta}.")
    else:
        seqrecord_printer(records)

# Execute if file ran as script
if __name__ == "__main__":
    main()


# Example of usage in Bash:
# python scr_random_dna_argparse.py -n 10 -l 100 -g 50 -i bashtest -f bashtest_randomdna.fasta