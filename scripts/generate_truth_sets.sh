#!/bin/bash

for i in tests/data/single_chrom; 
do
    python scripts/generate_truth_set.py -w 10 -o ${i} ${i}.fa ${i}.vcf.gz
done