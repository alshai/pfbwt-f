#!/bin/bash
if [[ -z ${1} ]]
then
    echo "no dir passed"
    exit 1
fi
SOURCE=${1}
mkdir -p out
python ${SOURCE}/vcf_to_bwt.py \
    --keep_parse \
    -o out/single_chrom \
    -m \
    -s \
    --ma_wsize 1 \
    --wsize 10 \
    ${SOURCE}/tests/data/single_chrom.fa \
    ${SOURCE}/tests/data/single_chrom.vcf.gz || { echo "vcf_to_bwt.py failed"; exit 1; }

python ${SOURCE}/scripts/readable_markers.py out/single_chrom.ma > out/single_chrom.markers
python ${SOURCE}/scripts/readable_sa.py out/single_chrom.sa > out/single_chrom.suffixarray
diff -qZ out/single_chrom.bwt ${SOURCE}/tests/data/single_chrom.bwt || { echo "BWT mismatch"; exit 1; }
diff -qZ out/single_chrom.suffixarray ${SOURCE}/tests/data/single_chrom.sa || { echo "SA mismatch"; exit 1; }
diff -qZ out/single_chrom.markers ${SOURCE}/tests/data/single_chrom.markers || { echo "marker mismatch exit 1"; exit 1; }
exit 0