#include <gtest/gtest.h>
#include "include/marker.hpp"


TEST(BitPacking, SimpleShift) {
    uint64_t x = 1;
    auto y = (x >> 1) << 1;
    EXPECT_EQ(y, 0);
}

TEST(BitPacking, SimpleShift7) {
    uint64_t x = 0xFuL << 60;
    uint64_t y = 0xF000000000000000uL;
    EXPECT_EQ(x, y);
}

TEST(BitPacking, SetAndGetPos1) {
    uint64_t x = set_pos(0, 100);
    auto y = get_pos(x);
    EXPECT_EQ(y, 100);
}

TEST(BitPacking, SetAndGetPos2) {
    uint64_t x = set_pos(0, 0x0010000000000000uL);
    auto y = get_pos(x);
    EXPECT_EQ(y, 0);
}

TEST(BitPacking, SetandGetSeq1) {
    uint64_t x = set_seq(0, 52);
    auto y = get_seq(x);
    EXPECT_EQ(y, 52);
}

TEST(BitPacking, SeqAndGetSeq2) {
    uint64_t x = set_seq(0, 0x10000);
    auto y = get_seq(x);
    EXPECT_EQ(y, 0);
}

TEST(BitPacking, SetandGetAle1) {
    uint64_t x = set_allele(0, 1);
    auto y = get_allele(x);
    EXPECT_EQ(y, 1);
}

TEST(BitPacking, SetAndGetAle2) {
    uint64_t x = set_allele(0, 0x10);
    auto y = get_allele(x);
    EXPECT_EQ(y, 0);
}

TEST(BitPacking, SetAndGetMarker1) {
    uint64_t x = 0;
    x = set_pos(x, 2839742);
    x = set_seq(x, 52);
    x = set_allele(x, 1);
    EXPECT_EQ(get_pos(x), 2839742);
    EXPECT_EQ(get_seq(x), 52);
    EXPECT_EQ(get_allele(x), 1);
}

TEST(BitPacking, SetAndGetMarker2) {
    uint64_t x = 0;
    x = set_pos(x, 2839742);
    x = set_seq(x, 52);
    x = set_allele(x, 0x10);
    EXPECT_EQ(get_pos(x), 2839742);
    EXPECT_EQ(get_seq(x), 52);
    EXPECT_EQ(get_allele(x), 0);
}

TEST(BitPacking, SetAndGetMarker3) {
    uint64_t x = 0;
    x = set_seq(x, 52);
    x = set_allele(x, 1);
    x = set_pos(x, 0x0010000000000000uL);
    EXPECT_EQ(get_pos(x), 0);
    EXPECT_EQ(get_seq(x), 52);
    EXPECT_EQ(get_allele(x), 1);
}