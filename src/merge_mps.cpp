#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>
#include "marker.hpp"
#include <map>

struct Args {
    std::string faidx;
    std::vector<std::string> fnames;
    std::string output;
    // uint64_t ref_length;
    int wsize = 10;
};

Args parse_args(int argc, char** argv) {
    Args args;
    if (argc < 5) {
        fprintf(stderr, "usage: ./merge_marker_indexes <faidx> <wsize> <output> <index 1> <index 2>");
        exit(1);
    }
    args.faidx = std::string(argv[1]);
    args.wsize = std::atoi(argv[2]);
    args.output.assign(argv[3]);
    for (int i = 4; i < argc; ++i) {
        args.fnames.push_back(std::string(argv[i]));
    }
    return args;
}


std::vector<size_t> lengths_from_faidx(std::string faidx) {
    FILE* faidx_fp = fopen(faidx.data(), "r");
    char* line = NULL;
    size_t n = 0;
    int nread;
    std::vector<size_t> lengths;
    while ((nread = getline(&line, &n, faidx_fp)) != -1) {
        std::string l(line);
        size_t pos = l.find("\t");
        std::string chr = l.substr(0, pos);
        l.erase(0, pos + 1);
        int len = std::atoi(l.substr(0, l.find("\t")).data());
        lengths.push_back(len);
    }
    return lengths;
}

void merge_indexes_mult(Args args) {
    auto lengths = lengths_from_faidx(args.faidx);
    size_t total_length = 0;
    for (auto l: lengths) total_length += l + args.wsize;
    FILE* ofp = fopen(args.output.data(), "wb");
    uint64_t x, mp, ms, pms=-1, ma;
    uint64_t delim = -1;
    int state = 0;
    size_t seq_bias = 0;
    uint64_t idx[2];
    std::vector<MarkerT> markers;
    int nseq = 0;
    for (auto f: args.fnames) fprintf(stderr, "%s\n", f.data());
    for (auto fname: args.fnames) {
        FILE* fp = fopen(fname.data(), "rb");
        if (fp == NULL) {
            fprintf(stderr, "error opening %s\n", fname.data());
            exit(1);
        }
        while (fread(&x, sizeof(uint64_t), 1, fp) == 1) {
            if (state < 2) { // keys are first two
                idx[state] = x;
                ++state;
            } else if (x != delim) { // values are rest
                ms = get_seq(x);
                markers.push_back(x);
                pms = ms;
            } else { // x == delim == separator between runs
                // deal with whole run here
                // fprintf(stderr, "%d %lu %lu %lu\n", nseq, ms, seq_bias, chrom_bias);
                idx[0] += seq_bias;
                idx[1] += seq_bias;
                fwrite(idx, sizeof(uint64_t), 2, ofp);
                fwrite(markers.data(), sizeof(MarkerT), markers.size(), ofp);
                fwrite(&x, sizeof(x), 1, ofp);
                markers.clear();
                state = 0;
            }
        }
        seq_bias += total_length;
        nseq += 1;
        pms = -1;
        fclose(fp);
    }
    fclose(ofp);
}

int main(int argc, char** argv) {
    merge_indexes_mult(parse_args(argc, argv));
    return 0;
}
