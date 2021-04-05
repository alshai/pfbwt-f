#ifndef MARKER_HPP
#define MARKER_HPP

#include <cinttypes>
#include <cassert>

using MarkerT = uint64_t;

constexpr uint64_t ALE_MASK = 0xF000000000000000uL;
constexpr uint64_t SEQ_MASK = 0x0FFFF00000000000uL;
constexpr uint64_t POS_MASK = 0x00000FFFFFFFFFFFuL;
constexpr int SEQ_SHIFT = 46;
constexpr int ALE_SHIFT = 60;

uint64_t set_pos(uint64_t x, uint64_t i) {
    return (x & ~POS_MASK) | (i & POS_MASK);
}

uint64_t get_pos(uint64_t x) {
    return (x & POS_MASK);
}

uint64_t set_seq(uint64_t x, uint64_t i) {
    return ((i & 0xFFFF) << SEQ_SHIFT) | (x & ~SEQ_MASK);
}

uint64_t get_seq(uint64_t x) {
    return (x & SEQ_MASK) >> SEQ_SHIFT;
}

uint64_t set_allele(uint64_t x, uint64_t i) {
    return ((i & 0xF) << ALE_SHIFT) | (x & ~ALE_MASK);
}

uint64_t get_allele(uint64_t x) {
    return (x & ALE_MASK) >> ALE_SHIFT;
}

inline MarkerT create_marker_t(uint64_t p, uint64_t a) {
    return set_pos(set_allele(0, a), p);
}

inline MarkerT create_marker_t(uint64_t pos, uint64_t ale, uint64_t seqid) {
    MarkerT x = 0;
    assert(ale < 16);
    assert(seqid < 65536); 
    assert(pos < 17592186044416);
    x = set_pos(x, pos);
    x = set_seq(x, seqid);
    x = set_allele(x, ale);
    return x;
}

// inline MarkerT create_marker_t(size_t pos, size_t ale, size_t seqid) {
//     MarkerT x = 0;
//     x = set_pos(x, static_cast<uint64_t>(pos));
//     x = set_seq(x, static_cast<uint64_t>(seqid));
//     x = set_allele(x, static_cast<uint64_t>(ale));
//     return x;
// }

#endif
