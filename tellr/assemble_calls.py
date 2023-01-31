import os
import gzip as gz

def de_novo(bam_name, repeat_fasta, sample, readfile):
    candidate_prefix = sample + '_candidates'
    print('samtools view -h {} -N {} -b -o {}.bam'.format(bam_name, str(readfile), candidate_prefix))
    os.system('samtools view -h {} -N {} -b -o {}.bam'.format(bam_name, readfile, candidate_prefix)) #save all candidate reads to bam

    print('samtools fasta  {}.bam > {}.fasta '.format(candidate_prefix, candidate_prefix))
    os.system('samtools fasta  {}.bam > {}.fasta '.format(candidate_prefix, candidate_prefix)) #bam to fasta

    paf_file = sample + '_candidate_overlap.paf.gz'
    print('minimap2 -x ava-pb {}.fasta {}.fasta| gzip -1 > {}'.format(candidate_prefix, candidate_prefix, paf_file))
    os.system('minimap2 -x ava-pb {}.fasta {}.fasta| gzip -1 > {}'.format(candidate_prefix, candidate_prefix, paf_file)) #map reads against themself, creating read overlap 

    gfa_file = sample + '_candidate_overlap.gfa'
    print('miniasm -f {}.fasta {} > {}'.format(candidate_prefix, paf_file, gfa_file))
    os.system('miniasm -f {}.fasta {} > {}'.format(candidate_prefix, paf_file, gfa_file)) #de novo assembly

    denovo_fasta = sample + '_denovo.fasta'
    gfa_fasta= """ awk '/^S/{{print ">"$2"\\n"$3}}' {0} > {1} """.format(gfa_file, denovo_fasta)
    os.system(gfa_fasta) #gfa to fasta

    repeat_paffile = sample + '_candidate_repeat_overlap.paf.gz'
    map_repeats = 'minimap2 -x ava-pb --cs {}  {} | gzip -1 >  {}'.format(denovo_fasta , repeat_fasta, repeat_paffile)
    print(map_repeats)
    os.system(map_repeats) #map de novo fasta to TE fasta
    

    return gfa_file, repeat_paffile


def breakpoints(gfa_file, repeat_paffile, readNamedict):
    
    #extracting de novo contig names and original read names
    denovoToReadName = {}
    for overlap in open(gfa_file):
        if not overlap.startswith('a') :
            continue

        denovo= overlap.split('\t')[1]
        reads = overlap.split('\t')[3].split(':')[0]
        if denovo not in denovoToReadName:
            denovoToReadName[denovo] = []
        denovoToReadName[denovo].append(reads)

        # paftools.js

    #repeats and contig matching more than 90% of repeat length
    #saved to dict denovo:repeat
    if 'gz' in repeat_paffile:
        opener = gz.open
    else:
        opener =open

    denovoToRepeat = {}
    for line in opener(repeat_paffile, 'rt'):

        r_len = int(line.split('\t')[1])
        match = int(line.split('\t')[9])
        perMatch = match/r_len

        if 'LINE' in line:
        #if perMatch > 0.7:
            repeat = line.split('\t')[0]
            denovo_contig = line.split('\t')[5]
            if denovo_contig not in denovoToRepeat:
                denovoToRepeat[denovo_contig] = []
            denovoToRepeat[denovo_contig].append(repeat)

    print(len(denovoToRepeat))
    for c in denovoToRepeat:
        print(c)
        reads = denovoToReadName[c]
        print(reads)
        for r in reads:
            pos = readNamedict[r]
            print(pos)

    return denovoToRepeat, denovoToReadName

def main(bam_name, repeat_fasta, sample, readfile, readName_dict):
    #aligned = de_novo(bam_name, repeat_fasta, sample, readfile)
    aligned = [sample + '_candidate_overlap.gfa', sample + '_candidate_repeat_overlap.paf.gz']
    breakpoints(aligned[0], aligned[1], readName_dict)
    #write vcf?

