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

def mergeDict(dict1, dict2):
    merged_dict = dict1.copy()
    merged_dict.update(dict2)
    return merged_dict
    
def main():
    parser = argparse.ArgumentParser(
        prog="mmf",
        description="Command-line script to perform motif finding of peaks.txt"
        
    )
    
    # Input
    parser.add_argument("peaks", help="HOMER peaks file", type=str)
    
    parser.add_argument("-f", "--fasta-ref", 
                        help="faidx Indexed Referencce Genome fasta file", 
                        metavar="FILE", type=str)
    
    parser.add_argument("-transfac", "--fac", help="transfac file from JASPAR", metavar="FILE", type=str)
    
    # Output
    parser.add_argument("-O", "--out", help="Write output to the directory." \
        "Default: stdout", metavar="DIR", type=str, required=False)
    
    # Parse args
    args = parser.parse_args()
    
    # genome.get_seq(start:, end:, chr)

    peaksDict = utils.getPeaksDict(args.peaks)
    pfm = utils.getFac(args.fac)
    genomeDict = loadGenome.load_genome(args.fasta_ref)
    # backgroundFreq = BF.getBackgroundFreq(genomeDict, peaksDict)
    # pwm = utils.getPWM(pfm)
    pwm = utils.getPWM(pfm, [0.28393081037720014, 0.2152101167620156, 0.2144817224377074, 0.2863773504230768])
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
    transpath = os.path.abspath(args.fac)
    
    if not os.path.exists(str(args.out)):
        os.makedirs(str(args.out))
    html = open(htmlpath, "w")
    
    header = """
<html>\n<head>\n<title> \nOutput Data in an HTML file
</title>\n</head> <h1>MMF Known Motif Enrichment Results </h1> \n
<h3> peaks.txt path: <u>{dir}</u> </h3>
<h2>TOP 10 Motif Found for <u>{transfac}</u></h2>\n
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
           border-collapse: collapse;}")

    html.write(header)
    n = 1
    for score in top10:
        html.write("\t<tr>\n")
        score_trimmed = format(score, '.5f')
        html.write("\t\t<th>{0}<th>{1}</th> <th>{2}</th>\n".format(n, scoresDict[score], score_trimmed))
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
