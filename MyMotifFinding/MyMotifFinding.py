#!/usr/bin/env python

"""
Command-line script to perform motif finding of peaks file

"""

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
    print(args.fasta_ref)
    print(args.out)
    print(args.peaks)
    print(args.fac)
    
    peaksDict = utils.getPeaksDict(args.peaks)
    pfm = utils.getFac(args.fac)
    genomeDict = loadGenome.load_genome(args.fasta_ref)
    # backgroundFreq = BF.getBackgroundFreq(genomeDict, peaksDict)
    pwm = utils.getPWM(pfm)
    # pwm = utils.getPWM(pfm, backgroundFreq)
    sequences = utils.getSequences(peaksDict, genomeDict)
    scoresDict = {}
    for seq in sequences:
        scoresDict = mergeDict(scoresDict, utils.getScore(pwm, seq))
    
    maxScore = max(list(scoresDict.keys()))
    bestMatchMotif = scoresDict[maxScore]
    print(bestMatchMotif)
    
    fasta = "C:\\Users\\Charles Choi\\Downloads\\GRCm38.chr17.fa"
    peaks = "C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\peaks.txt"
    fac = "C:\\Users\\Charles Choi\\Downloads\\MA0143.1.transfac"
if __name__ == "__main__":
    main()
