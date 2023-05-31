#!/usr/bin/env python

"""
Command-line script to perform motif finding of peaks file

"""
FASTA = "C:\\Users\\Charles Choi\\Downloads\\GRCm38.chr17.fa"
PEAKS = "C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\peaks.txt"
FAC = "C:\\Users\\Charles Choi\\Downloads\\MA0143.1.transfac"

import utils
import argparse
import os
import sys
import loadGenome
import Background_Frequency as BF
import findPValue as FP

def mergeDict(dict1, dict2):
    merged_dict = dict1.copy()
    merged_dict.update(dict2)
    return merged_dict
    
def main():
    parser = argparse.ArgumentParser(
        prog="mmf",
        description="Command-line script to perform motif finding of peaks.txt"
        
    )
    
    ## Input
    parser.add_argument("peaks", help="HOMER peaks file", type=str)
    
    parser.add_argument("-f", "--fasta-ref", 
                        help="faidx Indexed Referencce Genome fasta file", 
                        metavar="FILE", type=str)
    
    parser.add_argument("-transfac", "--fac", help="transfac file from JASPAR", metavar="FILE", type=str)
    
    ## Output
    parser.add_argument("-O", "--out", help="Write output to the directory." \
        "Default: stdout", metavar="DIR", type=str, required=False)
    
    ## Parse args
    args = parser.parse_args()
    
    ## Process data in peaks.txt
    peaksDict = utils.getPeaksDict(args.peaks)
    
    ## Process data in transfac file
    pfm = utils.getFac(args.fac)
    
    ## Process data in fasta reference genome
    genomeDict = loadGenome.load_genome(args.fasta_ref)
    
    ## Compute background frequency
    backgroundFreq = BF.getBackgroundFreq(genomeDict, peaksDict)
    
    ## Compute PWM from PFM
    pwm = utils.getPWM(pfm, backgroundFreq)
    
    ## Extract sequences from peaksDict
    sequences = utils.getSequences(peaksDict, genomeDict)
    scoresDict = {}
    for seq in sequences:
        scoresDict = mergeDict(scoresDict, utils.getScore(pwm, seq))
    
    scores = list(scoresDict.keys())
    scores.sort(reverse=True)
    top10 = scores[:10]
    
    ## html writing
    htmlpath = os.path.join(str(args.out), "KnownMotifFinding.html")
    peakDir = os.path.dirname(os.path.realpath(args.peaks))
    htmldir = os.path.abspath(args.out)
    transpath = os.path.abspath(args.fac)
    
    if not os.path.exists(str(args.out)):
        os.makedirs(str(args.out))
    html = open(htmlpath, "w")
    
    header = """
<html>\n<head>\n<title> \nOutput Data in an HTML file
</title>\n</head> <h1>MyMotifFinding Known Motif Enrichment Results </h1> <h2>(<u>{htmldir}</u>)</h2> \n
<h3> peaks.txt path: <u>{dir}</u> </h3>
<h3> transfac file path: <u>{transfac}</u></h3>\n
<style>
table, th, td {style}
</style>
<body>
<table style="width:100%">
    <tr>
        <th><b>Rank</b></th>
        <th><b>Motif</b></th>
        <th><b>Scores</b></th>
""".format(dir=peakDir, transfac=transpath, style="{border: 1px solid black;font-weight:400;\
           border-collapse: collapse;}", htmldir=htmldir)
    
    n = 1
    colorDict = {
        "A":"""<mark style="color:Red;background:none">""",
        "T":"""<mark style="color:Blue;background:none">""",
        "G":"""<mark style="color:Green;background:none">""",
        "C":"""<mark style="color:#ff8c00;background:none">"""
    }

    html.write(header)
    for score in top10:
        html.write("\t<tr>\n")
        score_trimmed = format(score, '.5f')
        html.write("\t\t<th>{0}</th><th><b>".format(n))
        for nuc in scoresDict[score]:
            html.write("{0}{1}</mark>".format(colorDict[nuc], nuc))
        html.write("</b></th><th>{0}</th>\n".format(score_trimmed))
        html.write("\t</tr>\n")
        n += 1
        
    tail = """
</table>
Thanks for using this tool!
</body></html>
    """
    
    html.write(tail)
    html.close()
    
    
if __name__ == "__main__":
    main()
run = """
python MyMotifFinding.py -f "C:\\Users\\Charles Choi\\Downloads\\GRCm38.chr17.fa" \
    -transfac "C:\\Users\\Charles Choi\\Downloads\\MA0143.1.transfac" \
        -O testFolder "C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\peaks.txt"
"""
