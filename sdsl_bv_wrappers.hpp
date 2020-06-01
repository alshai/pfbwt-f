#ifndef SDSL_BV
#define SDSL_BV

#include "sdsl/bit_vectors.hpp"

template<typename sdsl_bv_t=sdsl::bit_vector>
class bv_rs : public sdsl_bv_t {

    public:

    bv_rs() = default;

    bv_rs(const sdsl_bv_t& rhs) : sdsl_bv_t(rhs) { 
        init_rs();
    }

    bv_rs(sdsl_bv_t&& rhs) : sdsl_bv_t(rhs) { 
        init_rs();
    }

    bv_rs& operator=(const bv_rs& rhs) {
        sdsl_bv_t::operator=(rhs);
        init_rs();
        return *this;
    }

    bv_rs& operator=(bv_rs&& rhs) {
        sdsl_bv_t::operator=(rhs);
        init_rs();
        return *this;
    }

    bv_rs& operator=(const sdsl_bv_t& rhs) {
        sdsl_bv_t::operator=(rhs);
        init_rs();
        return *this;
    }

    bv_rs& operator=(sdsl_bv_t&& rhs) {
        sdsl_bv_t::operator=(rhs);
        init_rs();
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

    size_t size_in_bytes() {
        return sdsl::size_in_bytes(*this) + sdsl::size_in_bytes(rank_) + sdsl::size_in_bytes(select_);
    }

    private:

    typename sdsl_bv_t::rank_1_type   rank_;
    typename sdsl_bv_t::select_1_type select_;
};

#endif
