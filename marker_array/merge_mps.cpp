#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>

struct Args {
    std::vector<std::string> fnames;
    std::string output;
    uint64_t ref_length;
};

Args parse_args(int argc, char** argv) {
    Args args;
    if (argc < 5) {
        fprintf(stderr, "usage: ./merge_marker_indexes <ref length> <output> <index 1> <index 2>");
        exit(1);
    }
    args.ref_length = std::atol(argv[1]);
    args.output.assign(argv[2]);
    for (int i = 3; i < argc; ++i) {
        args.fnames.push_back(std::string(argv[i]));
    }
    return args;
}

void merge_indexes(Args args) {
    FILE* ofp = fopen(args.output.data(), "wb");
    int k = 0;
    // int w = 10;
    int64_t delta = 0;
    uint64_t pt = 0, pm = 0, x, delim = -1;
    int state = 0;
    for (auto fname: args.fnames) {
        FILE* fp = fopen(fname.data(), "rb");
        if (fp == NULL) {
            fprintf(stderr, "error opening %s\n", fname.data());
            exit(1);
        }
        while (fread(&x, sizeof(uint64_t), 1, fp) == 1) {
            if (state < 2) { // keys are first two
                uint64_t y = x + delta + (args.ref_length * k);
                fwrite(&y, sizeof(uint64_t), 1, ofp);
                pt = x;
                ++state;
            } else if (x != delim) { // values are rest
                fwrite(&x, sizeof(uint64_t), 1, ofp);
                pm = x;
            } else { // separator between runs
                fwrite(&delim, sizeof(uint64_t), 1, ofp);
                state = 0;
            }
        }
        delta += (pt - pm);
        // delta += w;
        ++k;
        fclose(fp);
    }
    fclose(ofp);
}

int main(int argc, char** argv) {
    merge_indexes(parse_args(argc, argv));
    return 0;
}
