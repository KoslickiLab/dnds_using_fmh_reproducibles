import pandas as pd
import argparse
from fmh_omega.helperfuncs import get_coding_sequence_from_nucleotide_sequence, positive_selection_based_on_mutation_rate_p, negative_selection_based_on_mutation_rate_p, read_fasta_gz
import random
from Bio.Seq import Seq
import gzip

def main(args):
    WD = Path(args.wd).expanduser().resolve()
    fasta_file = 'GB_GCA_021307345.1_protein.fna.gz
    mutation_rate_p = args.prate

    #get sequences from fasta file
    sequence_records = read_fasta_gz(fasta_file)

    # Create 100 random mutated sequences
    ITERATIONS = 100

    for i in range(ITERATIONS):
        #output files
        POSITIVE_SELECTION_QUERIES_NT = open(f'{WD}/positive_selection_queries_{mutation_rate_p}_{i}.fna','w')
        POSITIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}/positive_selection_translated_queries_{mutation_rate_p}_{i}.faa','w')
        NEGATIVE_SELECTION_QUERIES_NT = open(f'{WD}/negative_selection_queries_{mutation_rate_p}_{i}.fna','w')
        NEGATIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}/negative_selection_translated_queries_{mutation_rate_p}_{i}.faa','w')
        #loop through sequences of reference genome
        for name,ref_sequence in sequence_records.items():
            #ref sequence is mutated with mutation rate p
            query_positive_nt_seq = positive_selection_based_on_mutation_rate_p(ref_sequence,float(mutation_rate_p))
            POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_{mutation_rate_p}_{name}\n')
            POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

            #Translated mutated sequences 
            translated_positive_queries = Seq(query_positive_nt_seq).translate()
            POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_{mutation_rate_p}_{name}\n')
            POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_queries}\n')

            #ref sequence is mutated with mutation rate p
            query_negative_nt_seq = negative_selection_based_on_mutation_rate_p(ref_sequence,float(mutation_rate_p))
            NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_{mutation_rate_p}_{name}\n')
            NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

            #Translated mutated sequences 
            translated_negative_queries = Seq(query_negative_nt_seq).translate()
            NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_{mutation_rate_p}_{name}\n')
            NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_negative_queries}\n')

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate simulated protein coding sequence data based on positive and negative selection."
    )
    parser.add_argument(
        "--wd",
        type=str,
        default="./",
        help="Indicate working directory for output"
    )
    parser.add_argument(
        "--prate",
        type=float,
        default=0.01,
        help="Mutation rate p in [0,1]."
    )
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)