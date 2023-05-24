import scipy.stats
from Bio import motifs
from Bio.motifs.jaspar import calculate_pseudocounts
from Bio import SeqIO
import numpy as np

def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif is enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    peak_notb = peak_total-peak_motif
    bg_notb = bg_total - bg_motif
    table =[[peak_motif,bg_motif],[peak_notb,bg_notb]]
    odds, pval = scipy.stats.fisher_exact(table)
    
    return pval

def FindMaxScore(pwm, seq):
    """ Return the maximum score of a motif within a sequence """
    return max(pwm.calculate(seq))

# Load the PWMs from JASPAR
PWMList = [motifs.read(file, "jaspar") for file in pwm_files]
pwm_thresholds = [calculate_pseudocounts(pwm) for pwm in PWMList]

# Load the sequences 
peak_file = input("Enter the peak sequence file name: ")
bg_file = input("Enter the background sequence file name: ")
peak_seqs = [str(rec.seq) for rec in SeqIO.parse(peak_file, 'fasta')]
bg_seqs = [str(rec.seq) for rec in SeqIO.parse(bg_file, 'fasta')]


# Calculate the enrichment for each motif
for i in range(len(PWMList)):
    pwm = PWMList[i]
    thresh = pwm_thresholds[i]
    num_peak_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
    num_bg_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
    pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
    print(f"PWM: {pwm.name}, {num_peak_pass}/{len(peak_seqs)} peaks, {num_bg_pass}/{len(bg_seqs)} background; p-val: {pval}")
