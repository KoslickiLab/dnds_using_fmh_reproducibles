#!/usr/bin/env python3

import subprocess
import tempfile
import os

def load_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary mapping sequence IDs to their sequences.
    The sequence ID is assumed to be the first token (after '>') in each header.
    """
    sequences = {}
    with open(file_path, 'r') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq_lines)
                # Extract the sequence ID from the header (first word after '>')
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        # Save the last sequence in the file.
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq_lines)
    return sequences

def main():
    # File paths for input files
    orthologs_file = "total_orthologs_between_genomes.tab"
    
    # Updated mapping:
    # Query ID File 1 corresponds to this FASTA file:
    fasta_file_query1 = "scaffold_GCA_017434555.cds.fna"
    # Query ID File 2 corresponds to this FASTA file:
    fasta_file_query2 = "scaffold_GCA_017406775.cds.fna"
    
    # Create output directory for alignments
    output_dir = "cds_alignments"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load sequences from the FASTA files
    seqs_query1 = load_fasta(fasta_file_query1)
    seqs_query2 = load_fasta(fasta_file_query2)
    
    print("Loaded {} sequences from {}".format(len(seqs_query1), fasta_file_query1))
    print("Loaded {} sequences from {}".format(len(seqs_query2), fasta_file_query2))

    # After loading seqs_query1 and seqs_query2, collect IDs of shared orthologs
    shared1 = set()
    shared2 = set()
    
    # Open the orthologs file and process line by line.
    # It is assumed that the first line is a header containing column names.
    with open(orthologs_file, 'r') as f:
        header = f.readline().rstrip().split('\t')
        try:
            idx1 = header.index("Query ID File 1")
            idx2 = header.index("Query ID File 2")
        except ValueError:
            print("Error: Could not locate required header columns in {}".format(orthologs_file))
            return
        
        # Process each subsequent line for sequence pairs.
        for line_num, line in enumerate(f, start=2):
            parts = line.rstrip().split('\t')
            if len(parts) <= max(idx1, idx2):
                print("Skipping line {}: Not enough columns".format(line_num))
                continue
            
            # Use the appropriate FASTA file based on the query ID mapping:
            # Query ID File 1 comes from fasta_file_query1 and
            # Query ID File 2 comes from fasta_file_query2.
            id_query1 = parts[idx1]
            id_query2 = parts[idx2]
            
            seq1 = seqs_query1.get(id_query1)
            seq2 = seqs_query2.get(id_query2)
            
            if seq1 is None:
                print("Warning (line {}): Sequence '{}' not found in {}".format(line_num, id_query1, fasta_file_query1))
                continue
            if seq2 is None:
                print("Warning (line {}): Sequence '{}' not found in {}".format(line_num, id_query2, fasta_file_query2))
                continue

            # Record this pair as a shared ortholog
            shared1.add(id_query1)
            shared2.add(id_query2)
            

            """            
            # Write the two sequences to a temporary FASTA file for alignment.
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp_file:
                pair_file = tmp_file.name
                tmp_file.write(">{}\n".format(id_query1))
                tmp_file.write("{}\n".format(seq1))
                tmp_file.write(">{}\n".format(id_query2))
                tmp_file.write("{}\n".format(seq2))
            
            # Build a name for the alignment output file.
            # You might want to sanitize the IDs if they contain problematic characters.
            alignment_filename = "{}_{}.aln".format(id_query1, id_query2)
            output_file = os.path.join(output_dir, alignment_filename)
            
            print("Aligning sequences '{}' and '{}' (temporary file: {})".format(id_query1, id_query2, pair_file))
            
            # Perform the alignment with ClustalW using the temporary FASTA file as input.
            # The "-OUTFILE" option directs the output alignment to a specified file.
            return_code = subprocess.call([
                "clustalw",
                "-INFILE=" + pair_file,
                "-OUTPUT=FASTA",
                "-OUTFILE=" + output_file,
                "-QUIET"
            ])
            if return_code != 0:
                print("Error: clustalw returned non-zero exit code for pair '{}' and '{}'".format(id_query1, id_query2))
            else:
                print("Alignment saved to: {}".format(output_file))
            
            # Clean up the temporary file.
            os.remove(pair_file)
            """

    #I got dnds reported for these
    shared2 = {
    "JAFQXV010000014.1_48",
    "JAFQXV010000014.1_50",
    "JAFQXV010000014.1_51",
    "JAFQXV010000014.1_52",
    "JAFQXV010000014.1_53",
    "JAFQXV010000014.1_29",
    "JAFQXV010000014.1_56",
    "JAFQXV010000014.1_57",
    }

    shared1 = {
    "JAFRFB010000086.1_22",
    "JAFRFB010000086.1_24",
    "JAFRFB010000086.1_25",
    "JAFRFB010000086.1_26",
    "JAFRFB010000086.1_27",
    "JAFRFB010000086.1_3",
    "JAFRFB010000086.1_30",
    "JAFRFB010000086.1_31",
    }

    # Write out only the shared ortholog sequences into two new FASTA files
    orth1_out = os.path.join(output_dir, "GCA_017434555_shared_orthologs.fasta")
    orth2_out = os.path.join(output_dir, "GCA_017406775_shared_orthologs.fasta")
    with open(orth1_out, 'w') as f1, open(orth2_out, 'w') as f2:
        for q1 in sorted(shared1):
            f1.write(f">{q1}\n{seqs_query1[q1]}\n")
        for q2 in sorted(shared2):
            f2.write(f">{q2}\n{seqs_query2[q2]}\n")
    print(f"Wrote {len(shared1)} sequences to {orth1_out}")
    print(f"Wrote {len(shared2)} sequences to {orth2_out}")


if __name__ == "__main__":
    main()
