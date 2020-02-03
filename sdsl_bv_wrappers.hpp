#ifndef SDSL_BV
#define SDSL_BV

#include "sdsl/bit_vectors.hpp"

template<typename sdsl_bv_t=sdsl::bit_vector>
class bv_rs : public sdsl_bv_t {

    public:

    bv_rs() = default;

    bv_rs(const sdsl_bv_t& rhs) : sdsl_bv_t(rhs) { }

    bv_rs(sdsl_bv_t&& rhs) : sdsl_bv_t(rhs) { }

    bv_rs& operator=(const sdsl_bv_t& rhs) {
        sdsl_bv_t::operator=(rhs);
        return *this;
    }

    bv_rs& operator=(sdsl_bv_t&& rhs) {
        sdsl_bv_t::operator=(rhs);
        return *this;
    }

    void init_rs() {
        sdsl::util::init_support(rank_,   this);
        sdsl::util::init_support(select_, this);
    }

    size_t rank(size_t i) const {
        return rank_(i);
    }

    size_t select(size_t i) const {
        return select_(i);
    }

    private:

    sdsl::bit_vector::rank_1_type   rank_;
    sdsl::bit_vector::select_1_type select_;
};

#endif
