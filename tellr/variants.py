import sys
import pysam
import numpy as np


"""
fetch split reads
sort after split read size
get sequence
compare sequence to input file
"""

def check_cigar(flag, cigar, threshold):
    """
    Reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    threshold: needs to be at least threshold amount of soft clipped / inserted bases
    """
    basesToAdd={}

    num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold]) # For each pair (flag,bases) in cigar checks 0: flag and 1:bases to see if desired flag and >threshold

    if num_flags == 0:
        return False
    
    num_rounds =[] # Lenght of this to be same as length of num_flags

    add = 0     
    for c in cigar:
        if c[0] != flag:
            if c[0] != 1: # Dont want to add bases for insertions 
                add += c[1]

        else: # Matches flag
            if c in num_rounds: # Already counted
                add += c[1]
                continue

            elif int(c[1]) > threshold: # Threshold                 
                num_rounds.append(c)
            
                if add not in basesToAdd:
                    basesToAdd[add]=[]
                    basesToAdd[add].append(c[1])
                else:
                    basesToAdd[add].append(c[1])

                if len(num_rounds) == num_flags:
                    return basesToAdd
            else:
                if c[0]!= 1: # Dont want to add bases for insertions as output is reference base
                    add += c[1]
                continue

    return basesToAdd

def get_region(bamfile, chr, start, end, cand_list, cand_id_dict, cand_start_end, readstarts, mapping_quality, haploinfo, ReadToVarPos):
    """
    Runs through bam file and extracts positions of split reads and insertions in primary reads with mapping quality > user implied quality
    """
    appended = 0
    for read in bamfile.fetch(chr, start, end, until_eof=True):

        # Only check reads that are mapped and with good quality
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        elif read.mapping_quality < mapping_quality: # Mapping quality according to user set quality
            continue

        # Retrive HP tag for phasing/genotyping
        if read.qname not in haploinfo:
            haploinfo[read.qname] = 0
        try:
            haplotype = read.get_tag('HP') 
            haploinfo[read.qname] = haplotype
        except:
            haploinfo[read.qname] = 0

        cigg = read.cigar # Cigar: [ (1,24),(7,133) ]§
        read_start_pos = read.pos # Read Start

        # Keep track of which reads that are read
        if read.qname not in readstarts :
            readstarts[read.qname] = read_start_pos
        
        
        # Explain_cigar = 0: matching, 1:insertion, 2:deletion, 3:ref_skip, 4:soft_clipped, 5:hard_clipped, 6:...
        # Check for soft clipped 
        soft_clipping = check_cigar(4, cigg, 150) #list w bases to add for each soft clip OR false

        # Check for insertions
        insertions = check_cigar(1, cigg, 150) # list w bases to add for each ins OR false
 
        # Either could be a TE 
        if soft_clipping or insertions :

            # Find split read position
            # Reverse complement - opposite of others?
            
            variants = [soft_clipping, insertions]
            for v in variants:
                if v == False:
                    continue
                for add in v:
                
                    varpos = read_start_pos + add
                    varpos_end=varpos + max(v[add]) 

                    # Add to list for clustering
                    cand_list.append(varpos)

                    # Add to dict with chr:pos:[readname]
                    if varpos not in cand_id_dict:
                        cand_id_dict[varpos] = []
                    cand_id_dict[varpos].append(read.qname)

                    # Save Length of variant
                    if varpos not in cand_start_end:
                        cand_start_end[varpos] =[]
                    cand_start_end[varpos].append(varpos_end)

                    # Add name to dict with readname:pos
                    if read.qname not in ReadToVarPos:
                        ReadToVarPos[read.qname] = []
                    ReadToVarPos[read.qname].append(varpos)

                    appended += 1
            
        else:
            continue

    return cand_list, cand_id_dict, cand_start_end, readstarts, haploinfo, ReadToVarPos
 


def main(chr, bamfile, contig_length, mapping_quality):
    ReadToVarPos={}
    candidates = [] # [pos] return for each chromosome
    candidates_toid = {} # chr:{pos:[readname(s)]}
    cand_start_end = {} # Start and end of insert/split
    readstarts = {} # read:start
    haplotags = {} # readname : 1/2

    chr_len = contig_length[chr]
    for bin in range(1, chr_len-1000000, 1000000):
        get_region(bamfile, chr, bin, bin+1000000, candidates, candidates_toid , cand_start_end, readstarts, mapping_quality, haplotags,ReadToVarPos ) # Returns candidates (dict with split read and its supplementaries) for clustering

    return candidates, candidates_toid, cand_start_end, readstarts, haplotags, ReadToVarPos

