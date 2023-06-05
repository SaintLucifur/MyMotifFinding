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
import numpy as np
import time

FUNGI = os.path.join(os.getcwd(), "source", "fungi.transfac")
INSECTS = os.path.join(os.getcwd(), "source", "insects.transfac")
NEMATODES = os.path.join(os.getcwd(), "source", "nematodes.transfac")
PLANTS = os.path.join(os.getcwd(), "source", "plants.transfac")
UROCHORDATES = os.path.join(os.getcwd(), "source", "urochordates.transfac")
VERTEBRATES = os.path.join(os.getcwd(), "source", "vertebrates.transfac")

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
    
    parser.add_argument("-transfac", "--fac", help="taxonomic group from \
        fungi, insects, nematodes, plants, urochordates", type=str)
    
    ## Output
    parser.add_argument("-O", "--out", help="Write output to the directory." \
        "Default: stdout", metavar="DIR", type=str, required=False)
    
    ## Parse args
    args = parser.parse_args()
    
    ## Process data in peaks.txt
    peaksDict = utils.getPeaksDict(args.peaks)
    
    ## Get the correct transfac file
    transfacPathDict = {
        "fungi":FUNGI,
        "insects":INSECTS,
        "nematodes":NEMATODES,
        "plants":PLANTS,
        "urochordates":UROCHORDATES,
        "vertebrates":VERTEBRATES
    }
    TRANSFAC = transfacPathDict[args.fac]
    
    ## Process data in transfac file
    tf_time = time.time()
    id_pwm_logo_Dict = utils.getFac(TRANSFAC)
    print("\n--- transfac file processing: %.2f seconds ---\n" % (time.time() - tf_time))
    
    ## Process data in fasta reference genome
    gn_time = time.time()
    genomeDict = loadGenome.load_genome(args.fasta_ref)
    print("--- fasta file processing: %.2f seconds ---\n" % (time.time() - gn_time))
    
    ## Compute background frequency
    bg_time = time.time()
    backgroundFreq = BF.getBackgroundFreq(genomeDict, peaksDict)
    print("--- background freqs processing: %.2f seconds ---\n" % (time.time() - bg_time))
    
    ## Extract sequences from peaksDict
    sequences = utils.getSequences(peaksDict, genomeDict)
    
    ## Get PWM thresholds
    total = int(list(peaksDict.values())[0][4])
    if total*20 < 40000:
        numsim = total*20
    else:
        numsim = total
    
    ## Get number of PWMs
    pwmNum = len(id_pwm_logo_Dict.keys())
    print("--- total number of PWMs found: {0} ---\n".format(pwmNum))
    
    i = pwmNum
    
    pwm_time = time.time()
    for id in id_pwm_logo_Dict.keys():
        if i == pwmNum:
            first = time.time()
        bg_seqs = []
        pwm = id_pwm_logo_Dict[id][0]
        bg_seqs = [(utils.RandomSequence(pwm.shape[1], backgroundFreq)) for j in range(numsim)]
        null_scores = [utils.ScoreSeq(pwm, bg_seq) for bg_seq in bg_seqs]
        thresh = utils.GetThreshold(null_scores)
        num_peak_pass = np.sum([int(utils.FindMaxScore(pwm, seq)>thresh) for seq in sequences])
        num_bg_pass = np.sum([int(utils.FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
        pval = utils.ComputeEnrichment(total, num_peak_pass, numsim, num_bg_pass)
        id_pwm_logo_Dict[id].append("{:.0e}".format(pval))
        id_pwm_logo_Dict[id].append("{:.3e}".format(np.log10(-np.log10(pval+1)+1)))
        id_pwm_logo_Dict[id].append("{:.1f}".format(num_peak_pass))
        id_pwm_logo_Dict[id].append("{:.2f}".format(num_peak_pass/total*100))
        id_pwm_logo_Dict[id].append("{:.1f}".format(num_bg_pass))
        id_pwm_logo_Dict[id].append("{:.2f}".format(num_bg_pass/total*100))
        if i == pwmNum:
            totalTime = (time.time()-first)*pwmNum
            print("*** Estimated time: {:.2f} minutes ***\n".format(totalTime))
        i -= 1
        percentageDone = (1-i/pwmNum)*100

        if percentageDone != 100:
            print("--- percentage done {:.2f}% ---".format(percentageDone), end="\r")
        else:
            print("--- percentage done {:.2f}% ---\n".format(percentageDone))
        
    print("--- PWM processing: %.2f minutes ---\n" % ((time.time() - pwm_time)/60))
    print("*** Go ahead and copy the path of the output html file and open it in the browser! ***\n")
    
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
</title>\n</head> <h1>MyMotifFinding Known Motif Enrichment Results ({group})</h1> 
<h2>(<u>{htmldir}</u>)</h2> \n
<h3> peaks.txt path: <u>{dir}</u> </h3>
Total Target Sequences = {peaks}, Total Background Sequences = {numsim}
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
        <th>log p P-value</th>
        <th># peaks_match</th>
        <th>% peaks_match</th>
        <th># Bg_match</th>
        <th>% Bg_match</th>
""".format(dir=peakDir, style="{border: 1px solid black;font-weight:400;\
           border-collapse: collapse;}", htmldir=htmldir, numsim=numsim, peaks=total, group=args.fac)

    n = 1
    html.write(header)
    for tup in sortedTupleList:
        if n > 100:
            break
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


