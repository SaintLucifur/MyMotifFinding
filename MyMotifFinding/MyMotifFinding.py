#!/usr/bin/env python

"""
Command-line script to perform motif finding of peaks file

"""


import argparse
import os
import sys
# from pyfaidx import Fasta

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
    # Output
    parser.add_argument("-o", "--out", help="Write output to file." \
        "Default: stdout", metavar="FILE", type=str, required=False)
    
    # Parse args
    args = parser.parse_args()
    print(args.fasta_ref)
    print(args.out)
    print(args.peaks)
    
if __name__ == "__main__":
    main()