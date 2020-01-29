#!/bin/bash
INPUT=${1:-dengue.fa}
make pfbwt-f && \
rm -f tests/test/${INPUT}.* && \
./pfbwt-f -w 10 tests/test/${INPUT} && \
diff tests/test/${INPUT}.bwt tests/truth/${INPUT}.bwt && \
echo "success"
