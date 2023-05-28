import csv
import numpy as np
import math

fileName = ".\Test\peaksTest.txt"

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
        for row in f:
            if i < 39:
                i += 1
                continue
            else:
                peaksDict[row["#PeakID"]] = ["chr17", row["start"], row["end"],
                                                row["strand"]]
                
    peaksDict.pop("#PeakID", None)            
    return peaksDict

def getFac(facFile):
    pfmList = []
    with open(facFile, newline='') as csvFac:
        f = csv.reader(csvFac, delimiter='\t')
        i = 0
        for row in f:
            if i < 6:
                i += 1
                continue
            
            else:
                if row[0] == "XX":
                    break
                else:
                    pfmList.append(row)
                
    pfm = np.zeros((len(pfmList), 4))
    for i in range(len(pfm)):
        for j in range(len(pfm[0])):
            pfm[i][j] = pfmList[i][j+1]
    pfm = pfm.transpose()
    
    return pfm

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
        start = int(peaksDict[peak][1])
        end = int(peaksDict[peak][2])
        strand = peaksDict[peak][3]
        seq = genomeDict[chr][start:end]
        if strand == '-': ## handle the reverse read
            seq = getReverseComplement(seq)
        sequences.append(seq)
    return sequences

def getPWM(pfm, background_freqs=[0.25, 0.25, 0.25, 0.25]):

    pwm = np.zeros((4, len(pfm[0])))
    pfm = pfm + 0.01
    
    for i in range(len(pfm)):
        for j in range(len(pfm[0])):
            pwm[i][j] = math.log2(pfm[i,j]/np.sum(pfm[:,j])/background_freqs[i])
            
    return pwm

def getScore(pwm, seq):
    """return list of scores of motif on the given pwm using possible motifs from seq

    Args:
        pwm (numpy 2d array): positional weight matrix
        seq (str): long sequence of nucleotides

    Returns:
        float list: scores of motifs
    """
    n = pwm.shape[1]
    scoresDict = {}
    for i in range(len(seq)-n+1):
        scoresDict.update(ScoreSeq(pwm, seq[i:i+n]))
    return scoresDict

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
    scoreDict = {}
    
    nucs_rows = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(len(sequence)):
        score += pwm[nucs_rows[sequence[i]]][i]
        
    scoreDict[score] = sequence
    return scoreDict

# def main():
#     facFile = "C:\\Users\\Charles Choi\\Downloads\\MA0265.1.transfac"
#     pfm = getFac(facFile)
#     print(pfm)
    
#     pwm = getPWM(pfm)
#     print(pwm)
#     dict = getPeaksDict(fileName)
#     print(dict)
    
    ## Try Loading Genome
    # faFilePath = "C:\\Users\\Charles Choi\\Downloads\\hg38.fa"
    # genomeDict = LG.load_genome(faFilePath)
    # seqDict = getSequence(egPeaksDict, genomeDict)
    # print(seqDict)
    
