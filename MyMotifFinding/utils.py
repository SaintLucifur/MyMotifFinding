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


def getPFM(sequencesDict):
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pfm = np.zeros((4, len(list(sequencesDict.keys())[0])))
    
    for sequence in sequencesDict.keys():
        for j in range(len(sequence)):
            pfm[nucs[sequence[j]]][j] += sequencesDict[sequence]
        
    return pfm

# def getPFM(sequences):
#     nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
#     pfm = np.zeros((4, len(sequences[0])))
    
#     for i in range (len(sequences)):
#         for j in range(len(sequences[0])):
#             pfm[nucs[sequences[i][j]]][j] += 1
#     return pfm

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
    
    seqlogoCode = """
import seqlogo
# make seqlogo PWM object
seq_pwm = seqlogo.Pwm(pwm)
# Convert to ppm needed for plotting
seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
seqlogo.seqlogo(seq_ppm, ic_scale = True, format = 'png', size = 'medium')
    """
    
    jupy = nbformat.v4.new_notebook()
    pwmString = """pwm = """ + str(repr(np.array(pwm))) + '\n' + seqlogoCode
    print(pwmString)
    
    jupy['cells'] = [nbformat.v4.new_code_cell(pwmString)]
    nbformat.write(jupy,'Test.ipynb')
    
    # jupy['cells'] = [nbformat.v4.new_code_cell(seqlogoCode)]
    # nbformat.write(jupy,'Test.ipynb')
    
main()
