# TELLR

non-reference transposable element calling on long-read data 

# Run
require python, samtools, minimap2 and miniasm
python tellr/tellr.py [-h] [-ref REF] -TE_ref TE_REF -bam BAM


optional arguments:\n
  -h, --help            show this help message and exit \n
  -ref REF, --ref REF   reference genome \n
                        fasta sequence of elements to be detected \n
  -TE_ref TE_REF, --TE_ref TE_REF\n
  -bam BAM, --bam BAM   bam file\n