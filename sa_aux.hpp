#ifndef SA_AUX_HPP
#define SA_AUX_HPP

/* Author: Taher Mun
 * classes and functions for helping to build auxillary data structures to the
 * Suffix Array
 */

#include <cstdio>
#include <vector>
#include <tuple>
#include "sdsl/bit_vectors.hpp"
#include "sdsl_bv_wrappers.hpp"
extern "C" {
#include "utils.h"
}

namespace pfbwtf {

template<typename T>
class MarkerArray {

    public: 

    MarkerArray() {}

    MarkerArray(const char* fname) {
        std::FILE* fp = std::fopen(fname, "rb");
        std::fseek(fp, 0, SEEK_END);
        size_t fsize = std::ftell(fp);
        if (!fsize || (fsize % 3) || (fsize % sizeof(T))) {
            std::fprintf(stderr, "invalid marker file (size %lu)\n", fsize);
            std::exit(1);
        }
        std::fseek(fp, 0, SEEK_SET);
        size_t n = fsize / sizeof(T) / 3;

        std::vector<T> idx_v, allele_v;
        pos_v.reserve(n);
        idx_v.reserve(n);
        allele_v.reserve(n);
        T max_idx = 0;

        T x, y, z; // index, marker, allele
        while (fread(&x, sizeof(T), 1, fp) == 1 && // overall pos in text
            fread(&y, sizeof(T), 1, fp) == 1 && // pos w/i standard reference
            fread(&z, sizeof(T), 1, fp) == 1) { // allele idx (ref/alt)
            if (x > max_idx) { 
                max_idx = x;
            }
            idx_v.push_back(x);
            pos_v.push_back(y);
            if (z > 1) {
                die("sites w/ >1 alt alleles are not supported yet!");
            }
            allele_bv.push_back(z);
        }
        idx_bv.resize(max_idx+1);
        allele_bv.resize(max_idx+1);
        for (size_t i = 0; i < n; ++i) {
            idx_bv[idx_v[i]] = 1;
            allele_bv[idx_v[i]] = allele_v[i];
        }
        idx_bv.init_rs();
        allele_bv.init_rs();
    }

    std::pair<T, T> operator[](size_t i) const {
        if (idx_bv[i]) {
            return std::make_pair(pos_v[idx_bv.rank(i)], allele_bv[i]);
        }
        else return std::make_pair(0,0);
    }

    private: 

    bv_rs<> idx_bv;
    bv_rs<> allele_bv;
    std::vector<T> pos_v;
};

}; // end pfbwtf
#endif
