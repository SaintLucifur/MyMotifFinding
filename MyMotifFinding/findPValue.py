import numpy as np
from scipy.stats import fisher_exact
from Bio import SeqIO

def count_sequences(fasta_file, target_sequence):
    """Count the occurrence of target_sequence in the given fasta_file"""
    count = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        count += str(record.seq).count(target_sequence)
    return count

# Specify your target sequence and input files
target_sequence = "ATCG" # Replace this with your target sequence
file1 = "group1.fasta" # Replace with your actual file name
file2 = "group2.fasta" # Replace with your actual file name

# Count the target_sequence in each file
count1 = count_sequences(file1, target_sequence)
count2 = count_sequences(file2, target_sequence)

# Calculate the total number of sequences in each file
total1 = sum(1 for _ in SeqIO.parse(file1, "fasta"))
total2 = sum(1 for _ in SeqIO.parse(file2, "fasta"))

# Prepare a 2x2 contingency table
# table = [[count in file1, count in file2],
#          [total - count in file1, total - count in file2]]
table = np.array([[count1, count2], [total1 - count1, total2 - count2]])

# Perform Fisher's exact test
odds_ratio, p_value = fisher_exact(table)

# Print the results
print("Odds ratio:", odds_ratio)
print("p-value:", p_value)
