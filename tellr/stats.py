import pysam

"""
assess mapping quality, read depth -
in order to find regions with several split reads

"""


def quality(bamfile,chr, start, end ):
    for pileupcolumn in bamfile.pileup(chr, start , end):
        for pileupread in pileupcolumn.pileups:
            coverage = pileupcolumn.pos, pileupcolumn.n


            #paftools.js