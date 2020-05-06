#include <cstdio>
#include <cstdlib>
#include <cinttypes>

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: ./merge_marker_indexes <ref length> <output prefix> <index 1> <index 2>");
        exit(1);
    }
    size_t ref_length = std::atol(argv[1]);
    FILE* ofp = fopen(argv[2], "wb");
    int k = 0;
    int64_t delta = 0;
    uint64_t pt = 0, pm = 0, x, delim = -1; 
    int state = 0;
    for (int i = 3; i < argc; ++i) {
        FILE* fp = fopen(argv[i], "rb");
        if (fp == NULL) {
            fprintf(stderr, "error opening %s\n", argv[i]);
            exit(1);
        }
        while (fread(&x, sizeof(uint64_t), 1, fp) == 1) {
            if (x == delim) {
                fwrite(&delim, sizeof(uint64_t), 1, ofp);
                state = !state;
            } else if (state) { // value
                fwrite(&x, sizeof(uint64_t), 1, ofp);
                pm = x;
            } else { // key
                uint64_t y = x + delta + (ref_length * k);
                fwrite(&y, sizeof(uint64_t), 1, ofp);
                pt = x;
            }
        }
        delta += (pt - pm);
        ++k;
        fclose(fp);
    }
    fclose(ofp);
    return 0;
}
