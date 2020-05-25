#include <string>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include "marker_array/marker_index.hpp"
#include "rle_window_array.hpp"

void write_marker_array(std::string mai_fname, std::string sa_fname, std::string output = "") {
    FILE* sa_fp = sa_fname == "-" ? stdin : fopen(sa_fname.data(), "rb");
    FILE* ofp = fopen(output == "" ? "out" : output.data(), "wb");
    fprintf(stderr, "opening marker index from %s\n", mai_fname.data());
    MarkerIndex<> mai(mai_fname);
    constexpr uint64_t delim = -1;
    uint64_t s;
    uint64_t i = 0;
    std::vector<uint64_t> markers, pmarkers, locs;
    while (fread(&s, sizeof(uint64_t), 1, sa_fp) == 1) {
        markers.clear();
        markers = mai.get_markers(s);
        if (!vec_eq(markers, pmarkers)) {
            if (pmarkers.size()) {
                fwrite(&locs.front(), sizeof(uint64_t), 1, ofp);
                fwrite(&locs.back(), sizeof(uint64_t), 1, ofp);
                for (auto m: pmarkers) {
                    fwrite(&m, sizeof(uint64_t), 1, ofp);
                }
                fwrite(&delim, sizeof(uint64_t), 1, ofp);
            }
            locs.clear();
        }
        locs.push_back(i);
        pmarkers = markers;
        ++i;
    }
    if (pmarkers.size()) {
        fwrite(&locs.front(), sizeof(uint64_t), 1, ofp);
        fwrite(&locs.back(), sizeof(uint64_t), 1, ofp);
        for (auto m: pmarkers) {
            fwrite(&m, sizeof(uint64_t), 1, ofp);
        }
        fwrite(&delim, sizeof(uint64_t), 1, ofp);
    }
    fclose(sa_fp);
    fclose(ofp);
    return 0;
}


template<template<typename> typename ReadConType=VecFileSource>
class MarkerArray : public rle_window_arr<ReadConType> {

    public:

    MarkerIndex() {}
    MarkerIndex(std::string fname) : rle_window_arr<ReadConType>(fname) {}

    bool has_markers(uint64_t i) const {
        return this->has_entry(i);
    }

    std::vector<uint64_t> get_markers(uint64_t i) const {
        return this->at(i);
    }
};
