
import SeqIO

def load_fasta(file_path):
    
    with open(file_genome, "r") as handle:
        # Read the sequences using SeqIO.parse
        sequences = SeqIO.parse(handle, "fasta")
        # Create an empty dictionary to store chromosome sequences
        chr2seq = {}
        # Iterate over the sequences and populate the dictionary
        for record in sequences:
            chrom = record.id
            seq = str(record.seq).upper()
            chr2seq[chrom] = seq

    return chr2seq





# In[ ]:




