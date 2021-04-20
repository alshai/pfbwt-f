#!/bin/bash
if [[ -z ${1} ]]
then
    echo "no dir passed"
    exit 1
fi
SOURCE=${1}

if [[ -z ${2} ]]
then
    echo "no test passed"
    exit 1
fi

if [[ ! -f ${SOURCE}/tests/data/${2}.fa ]] || [[ ! -f ${SOURCE}/tests/data/${2}.vcf.gz ]];
then 
    echo "${SOURCE}/tests/data/${2}.fa or ${SOURCE}/tests/data/${2}.vcf.gz does not exist"
    exit 1
fi
TEST=${2}

mkdir -p out
python ${SOURCE}/vcf_to_bwt.py \
    --keep_parse \
    -o out/${TEST}.no_markers \
    -s \
    --wsize 10 \
    ${SOURCE}/tests/data/${TEST}.fa \
    ${SOURCE}/tests/data/${TEST}.vcf.gz || { echo "vcf_to_bwt.py failed"; exit 1; }

python ${SOURCE}/scripts/readable_sa.py out/${TEST}.no_markers.sa > out/${TEST}.no_markers.suffixarray
diff -qZ out/${TEST}.no_markers.bwt ${SOURCE}/tests/data/${TEST}.bwt || { echo "BWT mismatch"; exit 1; }
diff -qZ out/${TEST}.no_markers.suffixarray ${SOURCE}/tests/data/${TEST}.sa || { echo "SA mismatch"; exit 1; }
exit 0