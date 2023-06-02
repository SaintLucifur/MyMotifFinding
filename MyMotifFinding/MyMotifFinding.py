#!/usr/bin/env python

"""
Command-line script to perform motif finding of peaks file

"""
FASTA = "C:\\Users\\Charles Choi\\Downloads\\GRCm38.chr17.fa"
PEAKS = "C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\peaks.txt"
FAC = "C:\\Users\\Charles Choi\\Documents\\GitHub\\MyMotifFinding\\MyMotifFinding\\Test\\Two.transfac"
import utils
import argparse
import os
import sys
import loadGenome
import Background_Frequency as BF
import numpy as np


TRANSFAC = os.path.join(os.getcwd(), "source", "JASPAR2022_CORE_non-redundant_pfms_transfac.txt")

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
    
    ## Output
    parser.add_argument("-O", "--out", help="Write output to the directory." \
        "Default: stdout", metavar="DIR", type=str, required=False)
    
    ## Parse args
    args = parser.parse_args()
    
    ## Process data in peaks.txt
    peaksDict = utils.getPeaksDict(args.peaks)
    
    ## Process data in transfac file
    id_pwm_logo_Dict = utils.getFac(TRANSFAC)
    
    ## Process data in fasta reference genome
    genomeDict = loadGenome.load_genome(args.fasta_ref)
    
    ## Compute background frequency
    backgroundFreq = BF.getBackgroundFreq(genomeDict, peaksDict)
    
    ## Extract sequences from peaksDict
    sequences = utils.getSequences(peaksDict, genomeDict)
    
    ## Get PWM thresholds
    total = int(list(peaksDict.values())[0][4])
    numsim = total
    
    i = 1
    for id in id_pwm_logo_Dict.keys():
        bg_seqs = []
        pwm = id_pwm_logo_Dict[id][0]
        bg_seqs = [(utils.RandomSequence(pwm.shape[1], backgroundFreq)) for j in range(numsim)]
        null_scores = [utils.ScoreSeq(pwm, bg_seq) for bg_seq in bg_seqs]
        thresh = utils.GetThreshold(null_scores)
        num_peak_pass = np.sum([int(utils.FindMaxScore(pwm, seq)>thresh) for seq in sequences])
        num_bg_pass = np.sum([int(utils.FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
        pval = utils.ComputeEnrichment(total, num_peak_pass, numsim, num_bg_pass)
        id_pwm_logo_Dict[id].append("{:.2e}".format(pval))
        id_pwm_logo_Dict[id].append("{:.2e}".format(np.log10(pval+1)))
        id_pwm_logo_Dict[id].append("{:.1f}".format(num_peak_pass))
        id_pwm_logo_Dict[id].append("{:.2f}".format(num_peak_pass/total*100))
        id_pwm_logo_Dict[id].append("{:.1f}".format(num_bg_pass))
        id_pwm_logo_Dict[id].append("{:.2f}".format(num_bg_pass/total*100))
        # print("#{0}: {1}".format(i, pval))
    
    ## Rank p-value
    tupleList = []
    for key, value in id_pwm_logo_Dict.items():
        tup = (value[2], key)
        tupleList.append(tup)
    sortedTupleList = sorted(tupleList, key=lambda tup: float(tup[0]))

    ## html writing
    htmlpath = os.path.join(str(args.out), "KnownMotifFinding.html")
    peakDir = os.path.dirname(os.path.realpath(args.peaks))
    htmldir = os.path.abspath(str(args.out))
    
    if not os.path.exists(str(args.out)):
        os.makedirs(str(args.out))
    html = open(htmlpath, "w")
    
    header = """
<html>\n<head>\n<title> \nOutput Data in an HTML file
</title>\n</head> <h1>MyMotifFinding Known Motif Enrichment Results </h1> <h2>(<u>{htmldir}</u>)</h2> \n
<h3> peaks.txt path: <u>{dir}</u> </h3>
<style>
table, th, td {style}
</style>
<body>
<table style="width:100%">
    <tr>
        <th>Rank</th>
        <th>Motif</th>
        <th>Name</th>
        <th>P-value</th>
        <th>log P-value</th>
        <th># peaks_match</th>
        <th>% peaks_match</th>
        <th># Bg_match</th>
        <th>% Bg_match</th>
""".format(dir=peakDir, style="{border: 1px solid black;font-weight:400;\
           border-collapse: collapse;}", htmldir=htmldir)

    n = 1
    html.write(header)
    for tup in sortedTupleList:
        id = tup[1]
        html.write("\t<tr>\n")
        html.write("\t\t<th>{0}</th>".format(n)) ## Rank
        html.write("\t\t<th>{0}</th>".format(id_pwm_logo_Dict[id][1])) ## svg
        html.write("\t\t<th>{0}</th>".format(id)) ## Name
        html.write("\t\t<th>{0}</th>".format(id_pwm_logo_Dict[id][2])) ## p-value
        html.write("\t\t<th>{0}</th>".format(id_pwm_logo_Dict[id][3])) ## log p-value
        html.write("\t\t<th>{0}</th>".format(id_pwm_logo_Dict[id][4])) ## #peaks match
        html.write("\t\t<th>{0}%</th>".format(id_pwm_logo_Dict[id][5])) ## %peaks match
        html.write("\t\t<th>{0}</th>".format(id_pwm_logo_Dict[id][6])) ## #bg match
        html.write("\t\t<th>{0}%</th>".format(id_pwm_logo_Dict[id][7])) ## %bg match
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


