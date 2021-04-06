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
    -o out/mult_chroms \
    -m \
    -s \
    --ma_wsize 1 \
    ${SOURCE}/tests/data/mult_chroms.fa \
    ${SOURCE}/tests/data/mult_chroms.vcf.gz || { echo "vcf_to_bwt FAIL"; exit 1; }

python ${SOURCE}/scripts/readable_markers.py out/mult_chroms.ma > out/mult_chroms.markers
python ${SOURCE}/scripts/readable_sa.py out/mult_chroms.sa > out/mult_chroms.suffixarray
diff -qZ out/mult_chroms.bwt ${SOURCE}/tests/data/mult_chroms.bwt || { echo "BWT mismatch"; exit 1; }
diff -qZ out/mult_chroms.suffixarray ${SOURCE}/tests/data/mult_chroms.sa || { echo "SA mismatch"; exit 1; }
diff -qZ out/mult_chroms.markers ${SOURCE}/tests/data/mult_chroms.markers || { echo "marker mismatch exit 1"; exit 1; }
exit 0