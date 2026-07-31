import pandas as pd
import argparse
from simulate import get_coding_sequence_from_nucleotide_sequence, positive_selection_based_on_mutation_rate_p, negative_selection_based_on_mutation_rate_p
import random
from Bio.Seq import Seq
from pathlib import Path


def main(args):
    
    WD = Path(args.wd).expanduser().resolve()
    str_len = args.len #sequence length
    mutation_rate_p = args.prate

    #output files
    RANDOM_REF = open(f'{WD}/ref_{str_len}.fna','w')
    RANDOM_REF_PROT = open(f'{WD}/ref_translated_{str_len}.faa','w')

    POSITIVE_SELECTION_QUERIES_NT = open(f'{WD}/positive_selection_queries_{str_len}_{mutation_rate_p}.fna','w')
    POSITIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}/positive_selection_translated_queries_{str_len}_{mutation_rate_p}.faa','w')

    NEGATIVE_SELECTION_QUERIES_NT = open(f'{WD}/negative_selection_queries_{str_len}_{mutation_rate_p}.fna','w')
    NEGATIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}/negative_selection_translated_queries_{str_len}_{mutation_rate_p}.faa','w')

    # Create 100 random mutated sequences
    ITERATIONS = 100

    #Create a random X nt long sequence for simulation
    REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))

    #Get coding sequence and filter stop codons
    REF_coding_sequence = ''.join(get_coding_sequence_from_nucleotide_sequence(REF))

    #Save nt ref sequence to the following files
    RANDOM_REF.write(f'>ref_gene\n')
    RANDOM_REF.write(f'{REF_coding_sequence}\n')

    POSITIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
    POSITIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

    NEGATIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
    NEGATIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

    #save translated ref sequence to the following files
    ref_seq_translated = Seq(REF_coding_sequence).translate()
    RANDOM_REF_PROT.write(f'>ref_gene\n')
    RANDOM_REF_PROT.write(str(ref_seq_translated)+"\n")

    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>ref_gene\n')
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(str(ref_seq_translated)+"\n")

    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>ref_gene\n')
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(str(ref_seq_translated)+"\n")

    for i in range(ITERATIONS):
        #ref sequence is mutated with mutation rate p
        query_positive_nt_seq = positive_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
        POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_{mutation_rate_p}_{i}\n')
        POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

        #Translated mutated sequences 
        translated_positive_queries = Seq(query_positive_nt_seq).translate()
        POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_{mutation_rate_p}_{i}\n')
        POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_queries}\n')

        #ref sequence is mutated with mutation rate p
        query_negative_nt_seq = negative_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
        NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_{mutation_rate_p}_{i}\n')
        NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

        #Translated mutated sequences 
        translated_negative_queries = Seq(query_negative_nt_seq).translate()
        NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_{mutation_rate_p}_{i}\n')
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
        "--len",
        type=int,
        default=10002,
        help="Sequence length (positive integer)."
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