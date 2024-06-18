import argparse
from Bio import SeqIO

"""
This script can be used to convert genomic interval ranges of gaps within a chromosomal assembly into 
a bed file. It takes two arguments:
1. input fasta file witn "N" denoting gaps
2. the name of the chromosome contained within the input file.

It will output a bed file containing the length and position of each gap, sequentially.
"""


# Set up argument parsing
parser = argparse.ArgumentParser(description="Process a fasta file to identify sequence and gap windows.")
parser.add_argument("input_file", type=str, help="Path to the input fasta file.")
parser.add_argument("chromosome", type=str, help="Name of the chromosome contained in the input file.")

args = parser.parse_args()

# Assign arguments to variables
input_file = args.input_file
chromosome = args.chromosome

# Do not change these variables
output_file = chromosome + "_windows.bed"  # Output fasta can be found in "chromosome"_windows.bed
seq_position = 0
seq_list = []
seq_dict = {}
ATGC_count = 0  # Sequence window
N_count = 0  # Gap window

# Parse the fasta file
fasta_sequences = SeqIO.parse(input_file, "fasta")  # Iterator object

# Process the fasta sequences
for seq_record in fasta_sequences:
    for i in seq_record.seq:
        seq_list.append(i)  # Reads the sequence into a list of individual characters

for i in range(len(seq_list)):
    seq_position += 1
    if seq_list[i] == "N":  # Case if the character is an N
        if ATGC_count != 0:
            seq_dict[i] = [chromosome, seq_position - ATGC_count, seq_position - 1, "sequence", ATGC_count]  # Sequence window is over - add this section to the dictionary
            ATGC_count = 0  # Reset the window count
        N_count += 1  # Increment the gap window
    else:  # Case if the character is A, T, G, or C
        if N_count != 0:
            seq_dict[i] = [chromosome, seq_position - N_count, seq_position - 1, "gap", N_count]  # Gap window is over - add this section to the dictionary
            N_count = 0  # Reset the window count
        ATGC_count += 1  # Increment the sequence window

if ATGC_count != 0:  # Closing cases
    seq_dict[seq_position + 1] = [chromosome, seq_position - ATGC_count, seq_position - 1, "sequence", ATGC_count]  # Closes the last sequence and writes it into the file

if N_count != 0:  # Closing cases
    seq_dict[seq_position + 1] = [chromosome, seq_position - N_count, seq_position - 1, "gap", N_count]  # Closes the last gap and writes it into the file

# Writes the information to a bed file
with open(output_file, "w") as fasta_file:
    # Header for the columns in the file
    fasta_file.write("Chromosome\tStart\tStop\tType\tLength\n")
    # Write each key in the dictionary following the bed format
    for i in seq_dict:
        fasta_file.write('\t'.join(str(v) for v in seq_dict[i]) + "\n")
