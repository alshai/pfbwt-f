#include <string>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include "marker_array/marker_index.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: ./marker_index_to_array <marker index> <sa> [<output prefix>] \n");
    }
    FILE* sa_fp, *ofp;
    if (!strcmp("-", argv[2])) {
        fprintf(stderr, "reading sa from stdin\n");
        sa_fp = stdin;
    } else {
        fprintf(stderr, "opening %s as sa\n", argv[2]);
        sa_fp = fopen(argv[2], "rb");
    }
    if (argc > 3) ofp = fopen(argv[3], "wb");
    else ofp = fopen("out", "wb");

    std::string mai_fname = argv[1];
    fprintf(stderr, "opening marker index from %s\n", mai_fname.data());
    MarkerIndex<> mai(mai_fname);
    uint64_t key = 9411229;
    auto ms_results = mai.get_markers(key);
    std::cerr << "succinct index results for   " << key << ": ";
    for (auto i: ms_results) {
        std::cerr << i << " ";
    }
    std::cerr << "\n";
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

