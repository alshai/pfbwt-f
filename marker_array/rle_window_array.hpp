#ifndef RLE_WINDOW_HPP
#define RLE_WINDOW_HPP

#include <vector>
#include "sdsl_bv_wrappers.hpp"
#include "file_wrappers.hpp"

template<template<typename> typename ReadConType=VecFileSource>
class rle_window_arr {

    public:

    rle_window_arr() {}

    rle_window_arr(const bv_rs<>& run_starts, const bv_rs<>& run_ends, const std::vector<uint64_t>& arr)
        : run_starts_(run_starts)
        , run_ends_(run_ends)
        , arr_(arr) {}

    rle_window_arr(bv_rs<>&& run_starts, bv_rs<>&& run_ends, std::vector<uint64_t>&& arr)
        : run_starts_(run_starts)
        , run_ends_(run_ends)
        , arr_(arr) {}

    rle_window_arr(std::string fname) {
        ReadConType<uint64_t> in_arr(fname);
        uint64_t size = get_last_position(in_arr) + 2;
        run_starts_.resize(size);
        run_ends_.resize(size);
        uint64_t rs = 0, re=0;
        uint64_t keys[2];
        std::vector<uint64_t> values;
        int state = 0;
        for (size_t i = 0; i < in_arr.size(); ++i) {
            if (in_arr[i] == delim_) {
                if (rs == keys[0] || re == keys[1]) {
                    fprintf(stderr, "warning equal run start at %lu or run end at %lu\n", rs, re);
                    exit(1);
                }
                run_starts_[keys[0]] = 1;
                run_ends_[keys[1]] = 1;
                arr_.push_back(values);
                values.clear();
                state = 0;
            }
            else if (state < 2) {
                keys[state++] = in_arr[i];
            } else {
                values.push_back(in_arr[i]);
            }
        }
        run_starts_.init_rs();
        run_ends_.init_rs();
    }

    /*
    rle_window_arr& operator=(const rle_window_arr& rhs) {
        run_starts_ = rhs.run_starts_;
        run_ends_ = rhs.run_ends_;
        arr_ = rhs.arr_;
    }

    rle_window_arr& operator=(rle_window_arr&& rhs) {
        run_starts_ = std::move(rhs.run_starts_);
        run_ends_ = std::move(rhs.run_ends_);
        arr_ = std::move(rhs.arr_);
        run_starts_.init_rs();
        run_ends_.init_rs();
    }
    */

    bool has_entry(uint64_t i) const {
        return (run_starts_.rank(i+1) == run_ends_.rank(i) + 1);
    }

    std::vector<uint64_t> at(uint64_t i) const {
        uint64_t srank = i+1 > run_starts_.size() - 1 ? run_starts_.rank(run_starts_.size()-1) : run_starts_.rank(i+1);
        uint64_t erank = i > run_ends_.size()-1 ? run_ends_.rank(run_ends_.size()-1) : run_ends_.rank(i);
        if (srank != erank + 1) return std::vector<uint64_t>();
        else return arr_[srank-1];
    }

    private:

    uint64_t get_last_position(ReadConType<uint64_t>& in_arr) const {
        uint64_t i = in_arr.size() - 2; // in_arr[-1] is delim_
        while (in_arr[i] != delim_) {
            --i;
        }
        return in_arr[i+2]; //  delim keys[0] keys[1] val1 val2 ... delim EOF
    }


    bv_rs<> run_starts_;
    bv_rs<> run_ends_;
    std::vector<std::vector<uint64_t>> arr_;
    uint64_t delim_ = -1;
};

#endif
