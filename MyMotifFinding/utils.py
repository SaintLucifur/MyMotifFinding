import csv
import numpy as np
import math
import nbformat

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

def getSequence(peaksDict, genomeDict):
    sequenceDict = {}
    
    for peak in peaksDict.keys():
        chr = peaksDict[peak][0]
        start = peaksDict[peak][1]
        end = peaksDict[peak][2]
        count = peaksDict[peak][4]
        seq = genomeDict[chr][start:end]
        sequenceDict[seq] = count
        
    return sequenceDict


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


def main():
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
    pfm = getPFM(example)
    print(pfm)
    pwm = getPWM(pfm)
    print(pwm)
    
main()
