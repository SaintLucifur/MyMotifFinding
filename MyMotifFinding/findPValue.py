import scipy.stats
from Bio import motifs
from Bio.motifs.jaspar import calculate_pseudocounts
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
def calculate_pwm_thresholds(PWMList):
    """Calculate pseudocounts for each PWM in a list
    
    Parameters
    ----------
    PWMList : list
        A list of Position Weight Matrices (PWMs)

    Returns
    -------
    pwm_thresholds : list
        A list of pseudocounts for each PWM in the input list
    """
    return [calculate_pseudocounts(pwm) for pwm in PWMList]

# Load the sequences 
genome_dict = load_genome("GRCm38.chr17.fa")
peaksDict = getPeaksDict("peaks.txt")
peak_seqs = getSeqList(peaksDict, genome_dict)
bg_seqs = getBackgroundFreq(genome_dict, peaksDict, total_regions=10000)

# Calculate the enrichment for each motif
for i in range(len(PWMList)):
    pwm = PWMList[i]
    thresh = pwm_thresholds[i]
    num_peak_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
    num_bg_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
    pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
    print(f"PWM: {pwm.name}, {num_peak_pass}/{len(peak_seqs)} peaks, {num_bg_pass}/{len(bg_seqs)} background; p-val: {pval}")
