from typing import Dict

def load_genome(filename): 
    """ 
    This method load reference genome and convert it into a dictionary
    """
    genome_dict = {}
    chrom = None
    with open(filename, "r") as sequences:
        for line in sequences:
            line = line.strip()
            
            if line.startswith('>'):
                if chrom is not None:
                    genome_dict[chrom] = ''.join(genome_dict[chrom]).upper()
                chrom = line[1:]
                if chrom.startswith("chr"):
                    chrom = chrom[3:]
                genome_dict[chrom] = []
            else:
                genome_dict[chrom].append(line)

        if chrom is not None:
            genome_dict[chrom] = ''.join(genome_dict[chrom]).upper()

    return genome_dict  


#References: 
#https://stackoverflow.com/questions/29333077/reading-a-fasta-file-format-into-python-dictionary 
#https://stackoverflow.com/questions/22698807/parse-fasta-sequence-to-the-dictionary



    


