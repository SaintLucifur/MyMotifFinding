import csv
import numpy as np
import math
import nbformat
import loadGenome as LG

fileName = ".\Test\peaksTest.txt"

def getPeaksDict(fileName):
    """ Read the peaks.txt csv file

    Args:
        fileName (str): path to peaks.txt

    Returns:
        dict: key type: str   #PeakID
            value type: str list [chr, start, end, strand, Normalized Tag Count]
    """
    peaksDict = {}
    with open(fileName, newline='') as csvPeak:
        f = csv.DictReader(csvPeak, delimiter='\t', fieldnames=["#PeakID", "chr", "start", "end",
                        "strand", "Normalized Tag Count"])
        i = 0
        for row in f:
            if i < 39:
                i += 1
                continue
            else:
                print(row["#PeakID"], row["chr"], row["start"], 
                      row["end"], row["strand"], row["Normalized Tag Count"])
                
                peaksDict[row["#PeakID"]] = [row["chr"], row["start"], row["end"],
                                                row["strand"], row["Normalized Tag Count"]]
                
    peaksDict.pop("#PeakID", None)            
    return peaksDict

def getFac(facFile):
    pfmList = []
    with open(facFile, newline='') as csvFac:
        f = csv.reader(csvFac, delimiter='\t')
        i = 0
        hitXX = False
        for row in f:
            if i < 6:
                i += 1
                continue
            
            if hitXX == True:
                break
            
            if row[0] == "XX":
                hitXX = True
                
            else:
                pfmList.append(row)
                
    pfm = np.zeros((len(pfmList), 4))
    for i in range(len(pfm)):
        for j in range(len(pfm[0])):
            pfm[i][j] = pfmList[i][j+1]
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

def getSequence(peaksDict, genomeDict):
    sequenceDict = {}
    
    for peak in peaksDict.keys():
        chr = peaksDict[peak][0]
        start = int(peaksDict[peak][1])
        end = int(peaksDict[peak][2])
        count = int(peaksDict[peak][4])
        seq = genomeDict[chr][start:end]
        sequenceDict[seq] = count
        
    return sequenceDict

def getMotifDict(seqsDict, n):
    motifDict = {}
    
    for seq in seqsDict.keys():
        count = seqsDict[seq]
        for j in range(len(seq)-n+1):
            motif = seq[j:j+n]
            motifDict[motif] = count
            
    return motifDict

def getPFM(sequencesDict):
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pfm = np.zeros((4, len(list(sequencesDict.keys())[0])))
    
    for sequence in sequencesDict.keys():
        for j in range(len(sequence)):
            pfm[nucs[sequence[j]]][j] += sequencesDict[sequence]
        
    return pfm

def getBackgroundFreq():
    background_freq = [0.25, 0.25, 0.25, 0.25]
    return background_freq

def getPWM(pfm, background_freqs=[0.25, 0.25, 0.25, 0.25]):
    
    pwm = np.zeros((4, len(pfm[0])))
    pfm = pfm + 0.01
    
    for i in range(len(pfm)):
        for j in range(len(pfm[0])):
            pwm[i][j] = math.log2(pfm[i,j]/np.sum(pfm[:,j])/background_freqs[i])
            
    return pwm

def getScore(pwm, seq):
    n = pwm.shape[1]
    scores = [0]*(len(seq)-n+1) # list of scores. scores[i] should give the score of the substring sequence[i:i+n]
    # your code here
    for i in range(len(seq)-n+1):
        scores[i] = ScoreSeq(pwm, seq[i:i+n])
#     raise NotImplementedError
    return scores

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

def main():
    facFile = "C:\\Users\\Charles Choi\\Downloads\\MA0265.1.transfac"
    print(getFac(facFile))
    
    dict = getPeaksDict(fileName)
    for key in dict.keys():
        if key == '17-14':
            print(dict[key][1], dict[key][2])
    
    example = {
        "ATTAG":5,
        "TCGCG":6,
        "GGTGG":7,
        "TATAG":8,
        "TACGG":4,
        "AGCCG":3
    }
    # pfm = getPFM(example)
    # print(pfm)
    # pwm = getPWM(pfm)
    # print(pwm)
    motifDict = getMotifDict(example, 4)
    print(motifDict)
    
    egPeaksDict = {
        "Peak_1":["chr1", "1000", "1010", "+", "100"]
    }
    ## Try Loading Genome
    # faFilePath = "C:\\Users\\Charles Choi\\Downloads\\hg38.fa"
    # genomeDict = LG.load_genome(faFilePath)
    # seqDict = getSequence(egPeaksDict, genomeDict)
    # print(seqDict)
    
main()
