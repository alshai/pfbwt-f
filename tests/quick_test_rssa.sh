#!/bin/bash
INPUT=${1:-dengue.fa}
make pfbwt-f && \
rm -f tests/test/${INPUT}.* && \
./pfbwt-f -w 10 -r tests/test/${INPUT} && \
diff tests/test/${INPUT}.bwt tests/truth/${INPUT}.bwt && \
diff tests/test/${INPUT}.ssa tests/truth/${INPUT}.ssa && \
diff tests/test/${INPUT}.esa tests/truth/${INPUT}.esa && \
echo "success"
