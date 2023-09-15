import os
import gzip as gz
import pysam
import statistics

"""
maps sequence to TE and refines breakpoints
"""


def de_novo(chr, bam_name, repeat_fasta, sample, readfile, style):
    print('starting mapping')
    candidate_prefix = chr + '_' + sample + '_candidates'

    #extract reads of interest from bam to fasta
    bam_fasta = 'samtools view -h {} -N {} | samtools fasta -s stdout -n - > {}.fasta '.format(bam_name, str(readfile), candidate_prefix)
    print(bam_fasta)
    os.system(bam_fasta) 

    
    aligned_repeats = chr + '_' + sample + '_repeats.sam'

    #map fasta to TE fasta as short read 
    map_repeats = 'minimap2 -ax map-{} {}  {}.fasta >  {}'.format(style,repeat_fasta,candidate_prefix, aligned_repeats)
    print(map_repeats)
    os.system(map_repeats) 
    
    return aligned_repeats 



def check_cigar(flag, cigar, threshold):
    """
    reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    """
    basesToAdd=[]

    num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold])

    if num_flags == 0:
        return False
    num_rounds = []

    add = 0
    if cigar[0][0] == flag:
        add += int(cigar[0][1])
    for c in cigar: 
        if c[0] != flag:
            add += c[1]

        else:   
            if c[1] > threshold:
                if c in num_rounds:
                    add += c[1]
                    continue
                num_rounds.append(c)
                basesToAdd.append(add)
                if len(num_rounds) == num_flags:
                    return basesToAdd
            else:
                continue

    return basesToAdd

def breakpoints(chr, repeat_samfile,  sample, readNameToCluster, clusterToPos, readtopos ,PosToAvoid, readstarts):
    """
    finding breakpoint positions of TE
    """
    print("finding breakpoints")
    avoid_flags = [2048, 2064]

    repeat_vars = []
    #for line in open(repeat_samfile):
    samfile = pysam.AlignmentFile(repeat_samfile, 'r')
    samfile_header = samfile.header
    TEs = {}
    for te in samfile_header["SQ"]:
        thete = te["SN"]
        length = te["LN"]
        TEs[thete] = length
    
    for line in samfile.fetch():

        flag = line.flag
        if flag in avoid_flags:
            continue

        read = line.qname
        r_start = line.pos 
        repeat = line.reference_name
        if repeat not in TEs:
            continue
        #print(read, repeat)
        #connect start position in contig with position in read 
        cigar = line.cigar
        to_add = check_cigar(0, cigar, 100 )
        
        
        #find position in this read: read start + start of repeat

        cluster = list(readNameToCluster[read])
        for c in cluster:
            cluster_positions= list(clusterToPos[c])
            
            clusterconsensus = int(statistics.median(list(cluster_positions)))

            read_start = readstarts[read]
            
            #repeat_start = read_start + r_start
            repeat_start = readtopos[read]
            #print('repeat start', repeat_start)

            for r in repeat_start:
                pos_matches= [i for i in range(r-100, r+100) if i  in cluster_positions]
                if len(pos_matches) < 1:
                    #print('skip', read, repeat, r, pos_matches)
                    continue             
                    
                #check that its not a reference TE
                if chr in PosToAvoid:
                    start_positions = PosToAvoid[chr]
                    skip = [p for p in range(clusterconsensus-100, clusterconsensus+100) if p in start_positions ]#
                    if len(skip) > 0:
                        continue
                
                #print(c, cluster_positions, clusterconsensus)

                #print('read starts at', readstarts[read])
                #print('read split', readtopos[read])
                #print(read, flag, repeat, start, cluster, chr, clusterconsensus)
                lst = [repeat, chr, str(clusterconsensus), str(clusterconsensus + int(TEs[repeat])) ]
                #print(lst)
                if lst not in repeat_vars:
                    repeat_vars.append(lst)
    #print(repeat_vars)
    return repeat_vars


def pos_toavoid(file): #change to pytabix
    """
    poisitions to avoid (known reference TEs) provided by user
    """

    PosToAvoid = {}
    for line in open(file):
        if line.startswith('#'):
            continue
        chr = line.split('\t')[0].strip('chr')
        if chr not in PosToAvoid:
            PosToAvoid[chr] = {}
        apos = line.split('\t')[1]
        bpos = line.split('\t')[2]

        PosToAvoid[chr][apos] = bpos
    return PosToAvoid


def main(chr, bam_name, repeat_fasta, sample, readfile,  readNameToCluster, clusterToPos, read_topos, refrepeat, readstarts, style):
    aligned = de_novo(chr, bam_name, repeat_fasta, sample, readfile, style)
    #    aligned = chr + '_' + sample + '_repeats.sam' 
    avoid= pos_toavoid(refrepeat)
    #variants = breakpoints(chr, aligned[0], aligned[1], sample, readNameToCluster, clusterToPos, avoid)
    variants = breakpoints(chr, aligned, sample, readNameToCluster, clusterToPos, read_topos, avoid, readstarts)
  
    return variants

