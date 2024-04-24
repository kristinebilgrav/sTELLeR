import os
import gzip as gz
import pysam
import statistics

"""
maps sequence to TE and refines breakpoints
"""


def align(chr, bam_name, repeat_fasta, sample, readfile, style):
    print('starting mapping')
    
    aligned_repeats = chr + '_' + sample + '_repeats.sam'

    #map candidate-fasta to TE-fasta 
    map_repeats = 'minimap2 -ax map-{} {}  {} >  {}'.format(style,repeat_fasta,readfile, aligned_repeats)
    print(map_repeats)
    os.system(map_repeats) 
    
    return aligned_repeats 



def check_cigar(flag, cigar, threshold):
    """
    reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    """
    basesToAdd=0
    length=0
    match_start=False    
    interruptions=0
    num_flags = len([tup for tup in cigar if tup[0] == flag ])

    if num_flags == 0:
        return False
    num_rounds = []

    add = 0
    if cigar[0][0] == flag:
        add += int(cigar[0][1])
        length+=int(cigar[0][1])
    for c in cigar: 
        if c[0] != flag:

            if match_start:
                basesToAdd=add
                interruptions+=c[1]
                continue
            add += c[1]

        else:  
            match_start=True 
            length+=c[1]
            num_rounds.append(c)
            
            if num_rounds == num_flags:
                if length >= threshold : #and interruptions < (length/2)
                    return basesToAdd, length
            else:
                continue

    return basesToAdd, length


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


def extract_TEs(repeat_samfile):
    readToTEtype={} # Save readname and TE class
    readToTEPos={} # Save readname and position of identified TE

    avoid_flags=[4]

    samfile = pysam.AlignmentFile(repeat_samfile, 'r')
    samfile_header = samfile.header
    TEs = {}
    for te in samfile_header["SQ"]:
        thete = te["SN"]
        length = te["LN"]
        TEs[thete] = length
    
    for line in samfile.fetch():

        # Only check approved reads
        flag = line.flag
        if flag in avoid_flags:
            continue
        if line.mapping_quality < 20 or line.mapping_quality == 255 :
            continue
        if line.is_unmapped or line.is_duplicate:
            continue

        read = line.qname # Readname
        cigar=line.cigar
        match = check_cigar(0,cigar,100) # Returns where in read repeat starts and length of matching bases
        
        if match :
            repeat = line.reference_name # Repeat type
            if repeat not in TEs: # Confirm repeat in header
                continue

            # Save read to dictionaries
            if read not in readToTEtype:
                #readToTEPos[read]={}
                readToTEtype[read]=[] # Can be multiple TEs in one read
            readToTEtype[read].append(repeat) # Add all types of TEs found in read


            #if repeat not in readToTEPos[read]: # Save repeat start to repeat in read
               # readToTEPos[read][repeat] = [] # To verify the TE in same location as sr/ins
            #readToTEPos[read][repeat].append(tuple(match)) # Add bases and length


    return readToTEtype #, readToTEPos


def te_breakpoints(clusterToRead, readToTEtype,  chr, readinfo, sr):
    """
    Map confirmed TEs back to clusters
    Crosscheck that confirmed TE is same as clustered candidate
    //
    Go through each cluster 
    Check that reads in cluster are mapping to same TE
    Crosscheck that TE in read is same as clustered candidate
    """
  
    repeat_vars = []

    for c in clusterToRead: # go throguh each cluster and all reads associated with it
        tematches=[]
        reads = clusterToRead[c] # All reads in the cluster
        clusterinfo=[]

        # Go through reads in cluster and check their TE match and haplotag
        for r in reads:
            if r not in readToTEtype: # If read not in dict - did not match TE
                continue
           
            tes= readToTEtype[r] # Can have multiple TEs - list
            ht=readinfo[r][-2]
            VarStart= r.split('_')[-1] # Start of variant in read 
            varlen=readinfo[r][2]

            for te in tes:
                tematches.append(te)

            rinfo=[VarStart, ht, varlen]
            clusterinfo.append(rinfo)
           

        TEfrq=len(tematches)
        if TEfrq > sr:
            thete=max(set(tematches), key=tematches.count)

            pos= int(statistics.median([int(i[0]) for i in clusterinfo]))
            varlen=max([i[2] for i in clusterinfo])
            htags=set([str(i[1]) for i in clusterinfo])

            if len(htags) >1:
                teht=','.join(htags)
                tegt='1|1'
            else:
                teht=list(htags)[0]
                if teht == '1':
                    tegt='1|0'
                elif teht == '2':
                    tegt='0|1'
                else:
                    teht='0'
                    tegt='0/0'

                lst = [thete, chr, str(pos), str(pos + varlen), teht, tegt ]
                if lst not in repeat_vars:
                    repeat_vars.append(lst)

            # Check that its not a reference TE
            # if chr in PosToAvoid:
                # start_positions = PosToAvoid[chr]
                # skip = [p for p in range(clusterconsensus-100, clusterconsensus+100) if p in start_positions ]#
                # if len(skip) > 0:
                    # continue
 
    return repeat_vars


def main(chr, bam_name, repeat_fasta, sample, readfile, clusterToRead, refrepeat, style, readinfo, sr):
    aligned = align(chr, bam_name, repeat_fasta, sample, readfile, style)
    if os.path.isfile(refrepeat):
        avoid= pos_toavoid(refrepeat)

    tes=extract_TEs(aligned)
    readToTEtype =tes
    #readToTEPos = tes[1] 
    variants = te_breakpoints(clusterToRead, readToTEtype, chr, readinfo, sr)
    return variants

