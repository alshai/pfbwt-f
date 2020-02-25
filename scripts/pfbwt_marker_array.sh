#!/bin/bash

FA=${1:-"data/test.fa"};
VCF=${2:-"data/test.vcf.gz"};
# SAMPLES=${3:-"data/samples.txt"};
BYTES=8


# rm -f ${FA}.sa;
# mkfifo ${FA}.sa;

python scripts/vcf_haps_to_fasta.py --bytes 8 --stdout -o ${FA} ${FA} ${VCF}  |
    ./pfbwt-f64 -o ${FA} -s --stdout sa |
    python scripts/generate_marker_array.py -o ${FA} --bytes 8 - ${FA}.markers
