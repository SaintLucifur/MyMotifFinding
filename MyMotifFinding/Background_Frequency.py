#!/usr/bin/env python
# coding: utf-8

# In[1]:


from typing import Dict
import re
import random

def getSeqList(peaksDict, genome_dict):
    """
    return a list of sequences based on the peaks's start and end coordinates
    """
    sequenceList = []
    #chrom = returnChrom(peaksDict, genome_dict)
    for peak in peaksDict.keys():
        chr = peaksDict[peak][0]
        start = int(peaksDict[peak][1])
        end = int(peaksDict[peak][2])
        seq = genome_dict[chr][start:end]
        sequenceList.append(seq)
        
    return sequenceList


def find_peak_size_num(fileName):
    """
    Fidn the peak's size and number of peaks from the peaks.txt file
    """
    with open(fileName, 'r') as file:
        content = file.read()

    peakSize= None
    peakNum= None

    # Extract peak size
    peak_Match = re.search(r"# peak size = (\d+)", content)
    if peak_Match:
        peakSize= int(peak_Match.group(1))

    # Extract total peaks
    num_Match = re.search(r"# total peaks = (\d+)", content)
    if num_Match:
        peakNum= int(num_Match.group(1))

    return peakSize, peakNum

def calculate_GC(sequences):
    """
    Calculate GC content in a list of sequences
    """
    if len(sequences) == 0:
        return 0
    total_count = sequences.count("G")+sequences.count("C")+sequences.count("A")+sequences.count("T")
    if total_count == 0:
        return 0
    else:
        total_GC_count = sequences.count("G")+sequences.count("C")
        gc_percent = total_GC_count/total_count*100
    return gc_percent



def select_bg_regions(genome_dict, sequenceList, peakSize, peakNum, total_regions=10000):
    """
    Select background regions that match the GC distribution of the input list of sequences
    """
    load  = sequenceList
    chrom = returnChrom(peaksDict, genome_dict)
    load_seq = genome_dict[chrom]
    bg_regions_num = max(total_regions, 2 * peakNum)
    input_gc= [calculate_GC(seq) for seq in sequenceList if calculate_GC(seq) > 0]
    min_gc = min(input_gc)
    max_gc = max(input_gc)+5
    
    selected_bg_regions = []
    while len(selected_bg_regions) < bg_regions_num:
        start_pos = random.randint(0, len(load_seq)-peakSize)
        bg_region = load_seq[start_pos : start_pos + peakSize]
        bg_length = len(bg_region)
        bg_gcper = calculate_GC(bg_region)
        bg_N_per = bg_region.count("N")/bg_length * 100
        if bg_gcper >= min_gc and bg_gcper <= max_gc and bg_N_per < 70:
            selected_bg_regions.append(bg_region)
    return selected_bg_regions

    
def countFreq(selected_bg_regions):
    """
    Return a list containing the frequences of  ["A", "C", "G", "T"]
    """
    total_seq = "".join(selected_bg_regions)
    totalNum = total_seq.count("G")+total_seq.count("C")+total_seq.count("A")+total_seq.count("T")
    freqs = [0]*4
    freqs[0] = total_seq.count("A")/totalNum
    freqs[1] = total_seq.count("C")/totalNum
    freqs[2] = total_seq.count("G")/totalNum
    freqs[3] = total_seq.count("T")/totalNum
    return freqs


def getBackgroundFreq(genome_dict, peaksDict, total_regions=10000):
    """
    Return the background frequencies
    """
    peakNum = len(peaksDict) - 1 
    peakSize = None
    if peaksDict:
        fst_peak = next(iter((peaksDict.keys())))
        fst_v = peaksDict[fst_peak]
        peakSize = int(fst_v[2]) - int(fst_v[1]) 
    sequenceList = getSeqList(peaksDict, genome_dict)
    regions = select_bg_regions(genome_dict, sequenceList, peakSize, peakNum, total_regions)
    freq = countFreq(regions)
    return freq



# In[ ]:




