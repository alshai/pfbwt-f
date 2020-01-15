#!/bin/bash
make parse-f pfbwt-f && \
rm -f tests/test/dengue.fa.* && \
./parse-f -w 10 tests/test/dengue.fa && \
./pfbwt-f -w 10 tests/test/dengue.fa && \
diff tests/test/dengue.fa.bwt tests/truth/dengue.fa.bwt && \
echo "success"
