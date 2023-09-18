# TELLR

non-reference transposable element calling on long-read data 

# Run
require python3, pysam, samtools, minimap2 and miniasm

    python tellr [-h] [-R REF] -tf TE_FASTA -b BAM [-r SR] -tr TE_REF -s STYLE [-m MQ]

  arguments: 

    -R REF, --ref REF   reference genome

    -tf TE_FASTA, --TE_fasta TE_FASTA
                        fasta file with elements to be detected

    -b BAM, --bam BAM   bam file

    -tr TE_REF, --TE_ref TE_REF
                        bed file with positions to avoid

    -s STYLE, --style STYLE
                        ont or pb
  
  optional: 

    -sr SR, --sr SR       Minimum number of supporting split reads to call a variant (default 3)
    
    -maxsr

    -m MQ, --mq MQ        Mapping quality (default 20)
