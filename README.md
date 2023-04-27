# TELLR

non-reference transposable element calling on long-read data 

# Run
require python3, pysam, samtools, minimap2 and miniasm
  python tellr/tellr.py [-h] [-ref REF] -TE_fasta TE_FASTA -bam BAM -sr SR -TE_ref TE_REF

  arguments:
    -ref REF, --ref REF   reference genome
    -TE_fasta TE_FASTA, --TE_fasta TE_FASTA
                        fasta file with elements to be detected
    -bam BAM, --bam BAM   bam file
    -TE_ref TE_REF, --TE_ref TE_REF
                        bed file with positions to avoid
    -style STYLE, --style STYLE
                        ont or pb
  optional: 
    -sr SR, --sr SR       Minimum number of supporting split reads to call a variant (default 3)
    -maxsr