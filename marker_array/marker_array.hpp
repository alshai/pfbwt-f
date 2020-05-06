#include <cstdio>
#include <cinttypes>
#include <vector>
#include <string>
#include "parallel_hashmap/phmap.h"
#include "file_wrappers.hpp"

/* 
 * things I need
 * MarkerArrayWriter(idx_fname, sa_fname)
 *
 *
 */

struct MarkerArrayWriterArgs {
    std::string midx_fname;
    std::string sa_fname;
    int w = 10;
}

template< template<typename> typename ReadConType=VecFileSource>
class MarkerArrayWriter {

    using UIntType = uint64_t;

    public:
    
    using iterator = typename MarkerArr::iterator;

    MarkerArrayWriter() {}

    MarkerArrayWriter(MarkerArrayWriterArgs args) : w_(args.w) {
        MarkerIndex midx(args.midx_fname);
        FILE* sa_fp = fopen(args.sa_fname.data(), "rb");
        int key = 0;
        while (fread(&sa, sizeof(uint64_t), 1, sa_fp) == 1) {
            if (midx.find(sa) != midx.end()) {
                auto values = midx[sa];
            }
            ++key;
        }
        fclose(sa_fp);
    }

    void init_rs() {
        mbv.init_rs();
    }

    typename MarkerArr::const_iterator find(uint64_t i) {
        uint64_t j;
        if (mbv[i]) {
            j = mbv.rank(i);
        } else if (mbv.rank(i)) {
            j = mbv.rank(i)-1;
        } else return markers.cend();
        uint64_t midx = mbv.select(j+1);
        if (i >=midx && i < midx+markers[j].range_size) {
            return markers.cbegin() + j;
        } else {
            return markers.cend();
        }
    }

    // search a range for markers. return marker only if entire range contains it.
    typename MarkerArr::const_iterator find_range_single_only(uint64_t i1, uint64_t i2) {
        uint64_t j1, j2;
        if (mbv[i1]) {
            j1 = mbv.rank(i1);
        } else if (mbv.rank(i1)) {
            j1 = mbv.rank(i1)-1;
        } else return markers.cend();
        if (mbv[i2]) {
            j2 = mbv.rank(i2);
        } else if (mbv.rank(i2)) {
            j2 = mbv.rank(i2)-1;
        } else return markers.cend();
        if (j1 != j2) return markers.cend();

        uint64_t mstart = mbv.select(j1+1);
        uint64_t mend = mstart+markers[j1].range_size;
        if (i1 >= mstart && i1 < mend && i2 < mend) {
            return markers.cbegin() + j1;
        } else {
            return markers.cend();
        }
    }

    typename MarkerArr::iterator begin() {
        return markers.begin();
    }

    typename MarkerArr::const_iterator cbegin() const {
        return markers.begin();
    }

    typename MarkerArr::iterator end() {
        return markers.end();
    }

    typename MarkerArr::const_iterator cend() const {
        return markers.cend();
    }

    private:
    
    // TODO: turn this into one continuous vector
    MarkerArr markers;
    CompressedBv mbv;
};
}; // namespace marr

