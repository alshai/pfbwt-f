#!/bin/bash
INPUT=${1:-dengue.fa}
make parse-f pfbwt-f && \
rm -f tests/test/${INPUT}.* && \
./parse-f -w 10 -s tests/test/${INPUT} && \
./pfbwt-f -w 10 -s tests/test/${INPUT} && \
diff tests/test/${INPUT}.bwt tests/truth/${INPUT}.bwt && \
diff tests/test/${INPUT}.sa tests/truth/${INPUT}.sa && \
echo "success"
