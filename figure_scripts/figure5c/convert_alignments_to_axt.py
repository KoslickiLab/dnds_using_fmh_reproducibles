#!/usr/bin/env python3

import os
import subprocess

def fasta_to_sequences(fasta_file):
    """
    Reads a FASTA file and returns a list of tuples: (header, sequence).
    The sequence is concatenated into one continuous line.
    """
    sequences = []
    current_id = None
    current_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            sequences.append((current_id, "".join(current_seq)))
    return sequences

def main():
    # Directory where the FASTA alignment files are stored.
    alignments_dir = "cds_alignments"
    # Name of the combined AXT output file.
    combined_axt_file = "combined.axt"
    
    # Open the combined AXT file for writing.
    with open(combined_axt_file, 'w') as outfile:
        # Loop through each alignment file (assumed to have .aln extension).
        for filename in sorted(os.listdir(alignments_dir)):
            if not filename.endswith(".aln"):
                continue
            filepath = os.path.join(alignments_dir, filename)
            sequences = fasta_to_sequences(filepath)
            
            # Ensure that the file contains exactly two sequences.
            if len(sequences) != 2:
                print("Skipping file {}: expected 2 sequences but found {}"
                      .format(filename, len(sequences)))
                continue
            
            # Unpack the two sequences.
            header1, seq1 = sequences[0]
            header2, seq2 = sequences[1]
            
            # Create the combined AXT header using the two headers.
            combined_id = "{}_vs_{}".format(header1, header2)
            
            # Write the AXT entry: header line, first sequence, second sequence,
            # and an empty line as a separator.
            outfile.write("{}\n".format(combined_id))
            outfile.write("{}\n".format(seq1))
            outfile.write("{}\n".format(seq2))
            outfile.write("\n")
    
    print("Combined AXT file created: {}".format(combined_axt_file))
    
    # Define the KaKs_Calculator command.
    command = [
        "KaKs_Calculator", 
        "-i", combined_axt_file, 
        "-o", "kaks_results.txt", 
        "-m", "NG", 
        "-m", "LWL", 
        "-m", "LPB", 
        "-m", "GY", 
        "-m", "YN"
    ]
    
    print("Running KaKs_Calculator with command:")
    print(" ".join(command))
    
    # Run KaKs_Calculator using subprocess.
    return_code = subprocess.call(command)
    if return_code != 0:
        print("Error: KaKs_Calculator returned non-zero exit code {}".format(return_code))
    else:
        print("KaKs_Calculator successfully produced kaks_results.txt")

if __name__ == '__main__':
    main()
