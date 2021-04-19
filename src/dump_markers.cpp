#include <cstdio>
#include <cstdlib>
#include "marker.hpp"

int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: ./dump_markers <marker file>\n");
        exit(1);
    }
    FILE* fp = fopen(argv[1], "rb");
    if (fp == NULL) {
        fprintf(stderr, "error: file %s not found\n", argv[1]);
        exit(1);
    }
    uint64_t x;
    int state = 0;
    uint64_t delim = -1;
    while (fread(&x, sizeof(uint64_t), 1, fp)) {
        if (state < 2) {
            fprintf(stdout, "%lu\n", x);
            ++state;
        } else {
            if (x == delim) {
                fprintf(stdout, "%lu\n", x);
                state = 0;
            } else {
                fprintf(stdout, "%lu %lu %lu\n", get_seq(x), get_pos(x), get_allele(x));
            }
        }
    }
    return 0;
}
