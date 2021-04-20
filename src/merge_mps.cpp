#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>
#include "marker.hpp"
#include <map>

struct Args {
    std::vector<std::string> prefixes;
    std::string output;
};

Args parse_args(int argc, char** argv) {
    Args args;
    if (argc < 4) {
        fprintf(stderr, "usage: ./merge_marker_indexes <output> <prefix 1> ... <prefix n>");
        exit(1);
    }
    args.output.assign(argv[1]);
    for (int i = 2; i < argc; ++i) {
        args.prefixes.push_back(std::string(argv[i]));
    }
    return args;
}


void merge_indexes_mult(Args args) {
    FILE* ofp = fopen(args.output.data(), "wb");
    uint64_t x;
    uint64_t delim = -1;
    int state = 0;
    size_t seq_bias = 0;
    uint64_t idx[2];
    std::vector<MarkerT> markers;
    size_t length;
    for (auto prefix: args.prefixes) {
        std::string mps_fname = prefix + ".mps";
        std::string n_fname = prefix + ".n";
        FILE* mps_fp = fopen(mps_fname.data(), "rb");
        if (mps_fp == NULL) {
            fprintf(stderr, "error opening %s\n", mps_fname.data());
            exit(1);
        }
        /* get length of this sequence */
        FILE* n_fp = fopen(n_fname.data(), "r");
        if (n_fp == NULL) {
            fprintf(stderr, "error opening %s\n", n_fname.data());
            exit(1);
        } else {
            char* line = NULL;
            size_t n = 0;
            int nread;
            if ((nread = getline(&line, &n, n_fp)) != -1) {
                length = std::stoul(std::string(line));
            } else {
                fprintf(stderr, "error getting length from %s\n", n_fname.data());
                exit(1);
            }
            fclose(n_fp);
        }
        /* read in markers */
        while (fread(&x, sizeof(uint64_t), 1, mps_fp) == 1) {
            if (state < 2) { // keys are first two
                idx[state] = x;
                ++state;
            } else if (x != delim) { // values are rest
                markers.push_back(x);
            } else { // x == delim == separator between runs
                // deal with whole run here
                idx[0] += seq_bias;
                idx[1] += seq_bias;
                fwrite(idx, sizeof(uint64_t), 2, ofp);
                fwrite(markers.data(), sizeof(MarkerT), markers.size(), ofp);
                fwrite(&x, sizeof(x), 1, ofp);
                markers.clear();
                state = 0;
            }
        }
        seq_bias += length;
        fclose(mps_fp);
    }
    fclose(ofp);
}

int main(int argc, char** argv) {
    merge_indexes_mult(parse_args(argc, argv));
    return 0;
}
