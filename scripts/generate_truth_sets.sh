#!/bin/bash

for i in tests/data/single_chrom tests/data/mult_chroms tests/data/mult_chroms_indels; 
do
    python scripts/generate_truth_set.py -w 10 -o ${i} ${i}.fa ${i}.vcf.gz
done