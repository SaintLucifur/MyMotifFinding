#!/usr/bin/env python

"""
Command-line script to perform motif finding of peaks file

"""

from . import utils
import argparse
import os
import sys
import loadGenome

def mergeDict(dict1, dict2):
    return(dict2.update(dict1))

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
    pwm = utils.getPWM(pfm)
    sequences = utils.getSequences(peaksDict, genomeDict)
    scores_list = {}
    for seq in sequences:
        mergeDict(scores_list, utils.ScoreSeq(pwm, seq))
        
if __name__ == "__main__":
    main()