import scipy.stats
from Bio import motifs
from Bio.motifs.jaspar import calculate_pseudocounts
import numpy as np
from numpy.random import choice
import random

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

def RandomSequence(n, freqs):
    """ Generate a random string of nucleotides of length n
    
    Use the given nucleotide frequences
    
    Parameters
    ----------
    n : int
       Length of random string to generate
    freqs : list of float
       List of frequencies of A, C, G, T
       
    Returns
    -------
    seq : str
       random sequence of length n with the specified allele frequencies
    """
    seq = "A"*n
    seq = ""
    nucs = ['A', 'C', 'G', 'T']
    for i in range(n):
        seq += ''.join(random.choices(nucs, freqs))
    return seq

def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       pval% of null_dist should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 # set this  below to be the score threshold to obtain a p-value <0.01
    
    null_dist.sort(reverse=True)
    index = int(len(null_dist)*pval)
    thresh = null_dist[index]
    
    return thresh 
# def calculate_pwm_thresholds(PWMList):
#     """Calculate pseudocounts for each PWM in a list
    
#     Parameters
#     ----------
#     PWMList : list
#         A list of Position Weight Matrices (PWMs)

#     Returns
#     -------
#     pwm_thresholds : list
#         A list of pseudocounts for each PWM in the input list
#     """
#     return [calculate_pseudocounts(pwm) for pwm in PWMList]

# Load the sequences 


# def calculate_enrichment_for_each_pwm(PWMList, pwm_thresholds, peak_seqs, bg_seqs):
#     """Calculate the enrichment for each motif in the PWMList and print the result.

#     Parameters
#     ----------
#     PWMList : list
#         A list of Position Weight Matrices (PWMs)
#     pwm_thresholds : list
#         A list of pseudocounts for each PWM in the PWMList
#     peak_seqs : list
#         A list of peak sequences
#     bg_seqs : list
#         A list of background sequences
#     """
#     for i in range(len(PWMList)):
#         pwm = PWMList[i]
#         thresh = pwm_thresholds[i]
#         num_peak_pass = np.sum([int(FindMaxScore(pwm, seq) > thresh) for seq in peak_seqs])
#         num_bg_pass = np.sum([int(FindMaxScore(pwm, seq) > thresh) for seq in bg_seqs])
#         pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
#         print(f"PWM: {pwm.name}, {num_peak_pass}/{len(peak_seqs)} peaks, {num_bg_pass}/{len(bg_seqs)} background; p-val: {pval}")

# Call the function
# calculate_enrichment_for_each_pwm(PWMList, pwm_thresholds, peak_seqs, bg_seqs)

