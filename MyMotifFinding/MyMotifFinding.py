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
    parser.add_argument("-o", "--out", help="Write output to file." \
        "Default: stdout", metavar="FILE", type=str, required=False)
    
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
    html = open("KnownMotifFinding.html", "w")
    header = """
<html>\n<head>\n<title> \nOutput Data in an HTML file
</title>\n</head> <h1>MMF Known Motif Enrichment Results </h1> \n
<h3> ({dir}) </h3>
<h2>TOP 10 Motif Found for <u>dataset</u></h2>\n
<style>
table, th, td {style}
</style>
<body>
<table style="width:100%">
    <tr>
        <th>Motif</th>
        <th>Scores</th>
""".format(dir=os.getcwd(), style="{border:1px solid black;}")
    html.write(header)
    for score in top10:
        html.write("\t<tr>\n")
        html.write("\t\t<th>{0}</th> <th>{1}</th>\n".format(scoresDict[score], score))
        html.write("\t</tr>\n")
    tail = """
</table>
</body></html>
    """
    html.write(tail)
    html.close()
if __name__ == "__main__":
    main()
