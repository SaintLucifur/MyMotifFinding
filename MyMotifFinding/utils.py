import csv
import numpy as np
import math
from numpy.random import choice
import random
import scipy.stats

def getPeaksDict(fileName):
    """ Read the peaks.txt csv file

    Args:
        fileName (str): path to peaks.txt

    Returns:
        dict: key type: str   #PeakID
            value type: str list [chr, start, end, strand]
    """
    peaksDict = {}
    
    with open(fileName, newline='') as csvPeak:
        f = csv.DictReader(csvPeak, delimiter='\t', fieldnames=["#PeakID", "chr", "start", "end",
                        "strand"])
        i = 0
        for row in csvPeak:
            if i == 4:
                total = str(row).split(' ')[4]
                break
            i += 1
        
        i = 0
        for row in f:
            if i < 39:
                i += 1
                continue
            else:
                chr_value = row["chr"]
                if chr_value.startswith("chr"):
                    chr_value = chr_value[3:]
                peaksDict[row["#PeakID"]] = [chr_value, row["start"], row["end"],
                                                row["strand"], total]
                
    peaksDict.pop("#PeakID", None)            
    return peaksDict

def getFac(facFile, freqs=[0.25, 0.25, 0.25, 0.25]):
    """ parse the transfac file into PWMS and svgs

    Args:
        facFile (str): path to transfac file

        freqs (list): nucs frequencies
        
    Returns:
        dict: TF_ID : [PWM, svg]
    """
    pwmDict = {}
    with open(facFile, newline='') as csvFac:
        f = csv.reader(csvFac, delimiter='\t')
        i = 0
        for row in f:
            if row[0].startswith("DE"):
                pfmList = []
                id_asList = row[0].split(' ')
                de_ID = "{0} ID {1} {2}".format("from JASPAR:", id_asList[1], id_asList[2]) 
                row = next(f)
                row = next(f)
                while(row[0] != "XX"):
                    pfmList.append(row)
                    row = next(f)
                
                pfm = np.zeros((len(pfmList), 4))
                for i in range(len(pfm)):
                    for j in range(len(pfm[0])):
                        pfm[i][j] = pfmList[i][j+1]
                svg = getSVG(getPFM(pfm))
                pfm = pfm.transpose()
                
                pwm = getPWM(pfm, freqs)
                pwmDict[de_ID] = [pwm, svg]
                

    return pwmDict

def getReverseComplement(seq):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    
    reverse_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for nuc in seq[::-1]:
        revcomp += reverse_dict[nuc]
        
    return revcomp

def getSequences(peaksDict, genomeDict):
    sequences = []
    
    for peak in peaksDict.keys():
        chr = peaksDict[peak][0]
        if chr in genomeDict.keys():
            start = int(peaksDict[peak][1])
            end = int(peaksDict[peak][2])
            strand = peaksDict[peak][3]
            seq = genomeDict[chr][start:end]
            if strand == '-': ## handle the reverse read
                seq = getReverseComplement(seq)
            sequences.append(seq)
    return sequences

def getPFM(pfm):
    newpfm = np.zeros((len(pfm), 4))
    for i in range(len(pfm)):
        for j in range(4):
            newpfm[i][j] = pfm[i][j]/np.sum(pfm[i])
    return newpfm
            
def getPWM(pfm, background_freqs=[0.25, 0.25, 0.25, 0.25]):

    pwm = np.zeros((4, len(pfm[0])))
    pfm = pfm + 0.01
    
    for i in range(len(pfm)):
        for j in range(len(pfm[0])):
            pwm[i][j] = math.log2(pfm[i,j]/np.sum(pfm[:,j])/background_freqs[i])
            
    return pwm

###########################################################################################
## p-value calculation
def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    
    nucs_rows = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(len(sequence)):
        score += pwm[nucs_rows[sequence[i]]][i]
    
    return score

def FindMaxScore(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
            
    max_score = -1*np.inf
        
    n = len(pwm[0])
    scores = []
    for i in range(len(sequence)-n+1):
        scores.append(ScoreSeq(pwm, sequence[i:i+n]))
        scores.append(ScoreSeq(pwm, getReverseComplement(sequence[i:i+n])))
    if len(scores)>0:
        max_score = max(scores)
    return max_score

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

def GetThreshold(null_dist, pval=1e-5):
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

########################################################
## Get the image of motif logo
def getSVG(PFM):
    unitscript = """
<text fill="{0}" x="0" y="0"  transform="matrix({1},0,0,{2},{3},{4})">{5}</text>
    """
    script = """
<svg width="505" height="50">
    <g font-family="Arial" font-weight="bold" font-size="66.5">
    """
    nucDict = {
        0: "A",
        1: "C",
        2: "G",
        3: "T"
    }
    colorDict = {
        0: "#00BB00",
        1: "#0000EE",
        2: "#F9A500",
        3: "#DD0000"
    }
    fontDict = {
        0: 0.53,
        1: 0.59,
        2: 0.55,
        3: 0.65
    }
    plusDict = {
        0: 2,
        1: 0,
        2: 0,
        3: 1
    }
    
    n = 0
    for j in range(np.shape(PFM)[0]):
        posDict = {}
        compList = []
        for i in range(np.shape(PFM)[1]):
            propo = PFM[j][i]
            compList.append(propo)
            posDict[propo] = i
        
        sortedTuples = sorted(enumerate(compList), key=lambda i: i[1])
        assert(len(compList)==4)
        compList.sort()
        heightpos = 49
        for tuple in sortedTuples:
            index = tuple[0]
            propo = tuple[1]
            script += "\t" + unitscript.format(colorDict[index], fontDict[index], \
                propo, n+plusDict[index], heightpos, nucDict[index])
            heightpos = heightpos-50*propo
        n += 25
    script += "\t</g>\n</svg>"
    
    return script    
