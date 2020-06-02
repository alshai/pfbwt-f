#ifndef MARKER_HPP
#define MARKER_HPP

#include <cinttypes>

using MarkerT = uint64_t;

constexpr uint64_t ALE_MASK = 0xFF00000000000000;

constexpr uint64_t POS_MASK = 0x00FFFFFFFFFFFFFF;

inline uint8_t get_allele(MarkerT m) {
    return static_cast<uint8_t>(m >> 56);
}

inline uint64_t get_pos(MarkerT m) {
    return m & POS_MASK;
}

inline uint64_t set_allele(MarkerT m, uint8_t a) {
    return (m & POS_MASK) | (static_cast<uint64_t>(a) << 56);
}

inline uint64_t set_pos(MarkerT m, uint64_t p) {
    if ( p > POS_MASK ) {
        fprintf(stderr, "WARNING: %s: position %lu too large to pack into marker\n", __func__, p);
    }
    return (m & ALE_MASK) | p;
}

inline MarkerT create_marker_t(uint64_t p, uint64_t a) {
    return set_pos(set_allele(0, a), p);
}

#endif
