# TELLR

non-reference transposable element calling on long-read data 

# Run
require python3, pysam, samtools, minimap2 and miniasm
python tellr/tellr.py [-h] [-ref REF] -TE_fasta TE_FASTA -bam BAM -sr SR -TE_ref TE_REF


optional arguments:\n
  -h, --help            show this help message and exit \n
  -ref REF, --ref REF   reference genome \n
                        fasta sequence of elements to be detected \n
  -TE_ref TE_REF, --TE_ref TE_REF\n
  -TE_fasta TE_FASTA \n
  -bam BAM, --bam BAM   bam file\n