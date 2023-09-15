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
    reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    """
    basesToAdd=[]

    num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold])

    if num_flags == 0:
        return False
    num_rounds =[]
    
    

    add = 0
    if cigar[0][0] == flag:
        add += int(cigar[0][1])

     
    for c in cigar:
        if c[0] != flag:
            add += c[1]

        else: #matches flag
            if c in num_rounds: #already counted
                add += c[1]
                continue

            if c[1] > threshold: #threshold 
                
                num_rounds.append(c)
            
                basesToAdd.append(add)
                if len(num_rounds) == num_flags:
                    return basesToAdd
            else:
                continue

    return basesToAdd

def get_region(bamfile, samples, chr, start, end, cand_list, cand_id_dict, name_to_pos, readstarts, mapping_quality):
    """
    extracts split reads from bam file to dictionaries
    """
    appended = 0
    for read in bamfile.fetch(chr, start, end, until_eof=True):

        #only check reads that are mapped and with good quality
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        elif read.mapping_quality < mapping_quality:
            continue
        
        #get all split reads in region with support from more than X
        cigg = read.cigar
        read_start_pos = read.pos

        if read.qname not in readstarts :
            readstarts[read.qname] = read_start_pos
        

        #explain_cigar = 0: matching, 1:insertion, 2:deletion, 3:ref_skip, 4:soft_clipped, 5:hard_clipped, 6:...
        #check for soft clipped 
        soft_clipping = check_cigar(4, cigg, 150) #list w bases to add for each sfot clip OR false

        #check for insertions
        insertions = check_cigar(1, cigg, 150) # list w bases to add for each ins OR false
 
        #either could be a TE 
        if soft_clipping or insertions :

            #find split read position
            #reverse complement - oposite of others?
            
            variants = [soft_clipping, insertions]
            for v in variants:
                if v == False:
                    continue
                for add in v:
                
                    splitpos = read_start_pos + add

                    #if read.is_reverse:
                        #splitpos = read_start_pos - add

                #   add to list for clustering
                    cand_list.append(splitpos)

                    #add to dict with chr:pos:[readname]
                    if splitpos not in cand_id_dict:
                        cand_id_dict[splitpos] = []
                    cand_id_dict[splitpos].append(read.qname)

                    #add name to dict with readname:pos
                    if read.qname not in name_to_pos:
                        name_to_pos[read.qname] = []
                    name_to_pos[read.qname].append(splitpos)
                    appended += 1

            
        else:
            continue
    #print('total variants', appended)
    return cand_list, cand_id_dict, name_to_pos, readstarts
 



def main(chr, bamfile,bam_name, sample, contigs, contig_length, sr, mapping_quality):
    
    name_to_pos = {} #readname:chr-pos
    candidates = [] #[pos] return for each chromosome
    candidates_toid = {} #chr:{pos:[readname(s)]}
    readstarts = {} #read:start
    #if chr not in candidates:
            #candidates[chr] = []
    chr_len = contig_length[chr]
    for bin in range(1, chr_len-1000000, 1000000):
        get_region(bamfile, sample, chr, bin, bin+1000000, candidates, candidates_toid , name_to_pos, readstarts, mapping_quality) #returns candidates (dict with split read and its supplementaries) for clustering

    return candidates, candidates_toid, name_to_pos, readstarts

