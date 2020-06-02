#ifndef RLE_WINDOW_HPP
#define RLE_WINDOW_HPP

#include <vector>
#include "sdsl_bv_wrappers.hpp"
#include "file_wrappers.hpp"

template<template<typename> typename ReadConType=VecFileSource>
class rle_window_arr {

    public:

    rle_window_arr() {}

    /*
    rle_window_arr(const bv_rs<>& run_starts, const bv_rs<>& run_ends, const std::vector<std::vector<uint64_t>>& arr)
        : run_starts_(run_starts)
        , run_ends_(run_ends)
        , arr_(arr) {}

    rle_window_arr(bv_rs<>&& run_starts, bv_rs<>&& run_ends, std::vector<std::vector<uint64_t>>&& arr)
        : run_starts_(run_starts)
        , run_ends_(run_ends)
        , arr_(arr) {}
    */

    rle_window_arr(std::string fname) {
        ReadConType<uint64_t> in_arr(fname);
        uint64_t size = get_last_position(in_arr) + 2;
        sdsl::bit_vector run_starts, run_ends, arr_idxs;
        fprintf(stderr, "SIZE: %lu\n", size);
        run_starts.resize(size);
        run_ends.resize(size);
        arr_idxs.resize( in_arr.size() / (3 + wsize_) * wsize_ );
        uint64_t pos = 0;
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
                run_starts[keys[0]] = 1;
                run_ends[keys[1]] = 1;
                for (auto x: values) arr_.push_back(x);
                arr_idxs[pos] = 1;
                pos += values.size();
                // arr_.push_back(values);
                values.clear();
                state = 0;
            }
            else if (state < 2) {
                keys[state++] = in_arr[i];
            } else {
                values.push_back(in_arr[i]);
            }
        }
        run_starts_ = bv_t(run_starts);
        // run_starts_.init_rs();
        run_ends_ = bv_t(run_ends);
        // run_ends_.init_rs();
        arr_idxs.resize(arr_.size());
        arr_idxs_ = bv_t(arr_idxs);
        // arr_idxs_.init_rs();
    }

    void print_arr() {
        /*
        for (auto& v: arr_) {
            for (auto x: v) {
                fprintf(stderr, "%lu ", x);
            }
            fprintf(stderr, "\n");
        }
        */
        for (size_t i = 0; i < arr_.size(); ++i) {
            if (arr_idxs_[i] == 1) fprintf(stderr, "\n");
            fprintf(stderr, "%lu ", arr_[i]);
        }
        fprintf(stderr, "\n");
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

    /*
    std::vector<uint64_t> at(uint64_t i) const {
        uint64_t srank = i+1 > run_starts_.size() - 1 ? run_starts_.rank(run_starts_.size()-1) : run_starts_.rank(i+1);
        uint64_t erank = i > run_ends_.size()-1 ? run_ends_.rank(run_ends_.size()-1) : run_ends_.rank(i);
        if (srank != erank + 1) return std::vector<uint64_t>();
        else return arr_[srank-1];
    }
    */

    std::vector<uint64_t> at(uint64_t i, std::vector<uint64_t>& vals) const {
        uint64_t srank = i+1 > run_starts_.size() - 1 ? run_starts_.rank(run_starts_.size()-1) : run_starts_.rank(i+1);
        uint64_t erank = i > run_ends_.size()-1 ? run_ends_.rank(run_ends_.size()-1) : run_ends_.rank(i);
        if (srank != erank + 1) return std::vector<uint64_t>();
        fprintf(stderr, "srank: %lu - idxs: ", srank);
        uint64_t arr_idx_s = arr_idxs_.select(srank);
        fprintf(stderr, "%lu ", arr_idx_s);
        uint64_t arr_idx_e = srank + 1 > arr_idxs_.rank(arr_idxs_.size()) ? arr_idxs_.size() : arr_idxs_.select(srank + 1);
        fprintf(stderr, "%lu\n", arr_idx_e);
        vals.clear();
        for (size_t i = arr_idx_s; i < arr_idx_e; ++i) {
            vals.push_back(arr_[i]);
        }
        return vals;
    }

    std::vector<uint64_t> at(uint64_t i) const {
        std::vector<uint64_t> vals;
        return at(i, vals);
    }

    size_t calculate_arr_size() {
        arr_total_size_ = arr_.size();
        /*
        arr_total_size_ = 0;
        for (auto& v: arr_) {
            arr_total_size_ += (v.size() * sizeof(uint64_t));
        }
        return arr_total_size_;
        */
        return arr_total_size_;
    }

    size_t size_in_bytes() {
        if (!arr_total_size_) calculate_arr_size();
        // fprintf(stderr, "arr_: %lu\n", arr_total_size_);
        size_t run_starts_size = run_starts_.size_in_bytes();
        // fprintf(stderr, "run starts: %lu\n", run_starts_size);
        size_t run_ends_size = run_ends_.size_in_bytes();
        // fprintf(stderr, "run ends %lu\n", run_ends_size);
        return arr_total_size_ + run_starts_size + run_ends_size;
    }

    private:

    uint64_t get_last_position(ReadConType<uint64_t>& in_arr) const {
        uint64_t i = in_arr.size() - 2; // in_arr[-1] is delim_
        while (in_arr[i] != delim_) {
            --i;
        }
        fprintf(stderr, "last position at: %lu\n", i+2);
        return in_arr[i+2]; //  delim keys[0] keys[1] val1 val2 ... delim EOF
    }


    using bv_t = bv_rs<sdsl::sd_vector<>>;
    bv_t run_starts_;
    bv_t run_ends_;
    bv_t arr_idxs_;
    // std::vector<std::vector<uint64_t>> arr_;
    std::vector<uint64_t> arr_;
    uint64_t arr_total_size_ = 0;
    uint64_t delim_ = -1;
    int wsize_ = 10;
};

#endif
