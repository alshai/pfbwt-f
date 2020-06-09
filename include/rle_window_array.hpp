#ifndef RLE_WINDOW_HPP
#define RLE_WINDOW_HPP

#include <vector>
#include "sdsl_bv_wrappers.hpp"
#include "file_wrappers.hpp"

template<template<typename> typename ReadConType=VecFileSource>
class rle_window_arr {

    public:

    rle_window_arr() {}

    rle_window_arr(std::string fname) {
        ReadConType<uint64_t> in_arr(fname);
        uint64_t size = get_last_position_(in_arr) + 2;
        sdsl::bit_vector run_starts, run_ends, arr_idxs;
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
        run_ends_ = bv_t(run_ends);
        arr_idxs.resize(arr_.size());
        arr_idxs_ = bv_t(arr_idxs);
    }

    rle_window_arr(const rle_window_arr& rhs)
        : run_starts_(rhs.run_starts_)
        , run_ends_(rhs.run_ends_)
        , arr_idxs_(rhs.arr_idxs_)
        , arr_(rhs.arr_)
        , wsize_(rhs.wsize_)
    { }

    rle_window_arr(rle_window_arr&& rhs)
        : run_starts_(std::move(rhs.run_starts_))
        , run_ends_(std::move(rhs.run_ends_))
        , arr_idxs_(std::move(rhs.arr_idxs_))
        , arr_(std::move(rhs.arr_))
        , wsize_(rhs.wsize_)
    { }


    rle_window_arr& operator=(const rle_window_arr& rhs) {
        run_starts_ = rhs.run_starts_;
        run_ends_ = rhs.run_ends_;
        arr_idxs_ = rhs.arr_idxs_;
        arr_ = rhs.arr_;
        wsize_ = rhs.wsize_;
    }

    rle_window_arr& operator=(rle_window_arr&& rhs) {
        run_starts_ = std::move(rhs.run_starts_);
        run_ends_ = std::move(rhs.run_ends_);
        arr_idxs_ = std::move(rhs.arr_idxs_);
        arr_ = std::move(rhs.arr_);
        wsize_ = rhs.wsize_;
    }

    bool operator==(const rle_window_arr& rhs) const {
        if (run_starts_ != rhs.run_starts_) {
            return false;
        }
        if (run_ends_ != rhs.run_ends_) {
            return false;
        }
        if (arr_idxs_ != rhs.arr_idxs_) {
            return false;
        }
        if (arr_.size() != rhs.arr_.size()) {
            return false;
        }
        for (size_t i = 0; i < arr_.size(); ++i) {
            if (arr_[i] != rhs.arr_[i]) {
                return false;
            }
        }
        return true;
    }

    bool has_entry(uint64_t i) const {
        return (run_starts_.rank(i+1) == run_ends_.rank(i) + 1);
    }

    /* TODO: should we clear vals beforehand? */
    std::vector<uint64_t>& at(uint64_t i, std::vector<uint64_t>& vals) const {
        uint64_t srank = run_starts_rank(i); // this is the run-start before i
        uint64_t erank = run_ends_rank(i); // this is run-end before i
        // make sure that erank occurs before srank
        if (srank != erank + 1) return vals;
        return arr_at_(srank-1, vals);
    }

    std::vector<uint64_t> at(uint64_t i) const {
        std::vector<uint64_t> vals;
        return at(i, vals);
    }

    // find all markers within range [s,e], inclusive
    // TODO: reduce amount of select queries.
    // TODO: clear vals here?
    std::vector<uint64_t>& at_range(uint64_t s, uint64_t e, std::vector<uint64_t>& vals) const {
        uint64_t e_rs_rank = run_starts_rank(e);
        uint64_t e_rs_pos = run_starts_select(e_rs_rank);
        if (e_rs_pos <= s) {
            uint64_t e_re_pos = run_ends_select(e_rs_rank);
            return e_re_pos >= s ? arr_at_(e_rs_rank, vals) : vals;
        }
        uint64_t s_rs_rank = run_starts_rank(s);
        uint64_t s_rs_pos = run_starts_select(s_rs_rank);
        uint64_t s_re_rank = run_ends_rank(s);
        uint64_t s_re_pos = run_ends_select(s_re_rank);
        uint64_t start_idx = s_rs_pos > s_re_pos ? s_rs_rank - 1 : s_rs_rank;
        assert(run_starts_select(start_idx+1) <= e); // this check should have been taken care of above
        for (uint64_t i = start_idx; i < e_rs_rank; ++i) {
            arr_at_(i, vals);
        }
        return vals;
    }

    std::vector<uint64_t> at_range(uint64_t s, uint64_t e) const {
        std::vector<uint64_t> vals;
        return at_range(s, e, vals);
    }

    void print_arr() const {
        for (size_t i = 0; i < arr_.size(); ++i) {
            if (arr_idxs_[i] == 1) fprintf(stderr, "\n");
            fprintf(stderr, "%lu ", arr_[i]);
        }
        fprintf(stderr, "\n");
    }

    size_t serialize(std::ostream& out) {
        size_t nbytes = 0;
        nbytes += run_starts_.serialize(out);
        nbytes += run_ends_.serialize(out);
        nbytes += arr_idxs_.serialize(out);
        size_t arr_size = arr_.size();
        out.write(reinterpret_cast<char*>(&arr_size), sizeof(arr_size));
        nbytes += sizeof(arr_size);
        out.write(reinterpret_cast<char*>(arr_.data()), sizeof(arr_[0]) * arr_.size());
        nbytes += sizeof(arr_[0]) * arr_.size();
        out.write(reinterpret_cast<char*>(&wsize_), sizeof(wsize_));
        nbytes += sizeof(wsize_);
        return nbytes;
    }

    void load(std::istream& in) {
        run_starts_.load(in);
        run_ends_.load(in);
        arr_idxs_.load(in);
        size_t arr_size;
        in.read(reinterpret_cast<char*>(&arr_size), sizeof(arr_size));
        arr_.resize(arr_size);
        in.read(reinterpret_cast<char*>(arr_.data()), sizeof(uint64_t) * arr_size);
        in.read(reinterpret_cast<char*>(&wsize_), sizeof(wsize_));
    }

    private:

    inline uint64_t run_starts_rank(uint64_t i) const {
        return i+1 >= run_starts_.size()
                     ? run_starts_.rank(run_starts_.size())
                     : run_starts_.rank(i+1);
    }

    // will fail on i == 0
    inline uint64_t run_starts_select(uint64_t i) const {
        return i > run_starts_.rank(run_starts_.size())
               ? run_starts_.size()
               : run_starts_.select(i);
    }

    inline uint64_t run_ends_rank(uint64_t i) const {
         return i > run_ends_.size()-1
             ? run_ends_.rank(run_ends_.size()-1)
             : run_ends_.rank(i);
    }

    // will fail on i == 0
    inline uint64_t run_ends_select(uint64_t i) const {
        return i > run_ends_.rank(run_ends_.size())
               ? run_ends_.size()
               : run_ends_.select(i);
    }

    inline uint64_t arr_idxs_select(uint64_t i) const {
        return i > arr_idxs_.rank(arr_idxs_.size())
               ? arr_idxs_.size()
               : arr_idxs_.select(i);
    }

    // i is 0-index into theoretical arr of arrs
    // remember that rank(i) starts from 1, so always pass rank(i)-1
    std::vector<uint64_t>& arr_at_(uint64_t i, std::vector<uint64_t>& vals) const {
        uint64_t arr_s_idx = arr_idxs_select(i+1);
        uint64_t arr_e_idx = arr_idxs_select(i+2);
        for (uint64_t i = arr_s_idx; i < arr_e_idx; ++i) {
            vals.push_back(arr_[i]);
        }
        return vals;
    }

    std::vector<uint64_t> arr_at_(uint64_t i) const {
        std::vector<uint64_t> vals;
        return arr_at_(i, vals);
    }

    uint64_t get_last_position_(ReadConType<uint64_t>& in_arr) const {
        uint64_t i = in_arr.size() - 2; // in_arr[-1] is delim_
        while (in_arr[i] != delim_) {
            --i;
        }
        return in_arr[i+2]; //  delim keys[0] keys[1] val1 val2 ... delim EOF
    }

    using bv_t = bv_rs<sdsl::sd_vector<>>;
    bv_t run_starts_;
    bv_t run_ends_;
    bv_t arr_idxs_;
    std::vector<uint64_t> arr_;
    uint64_t delim_ = -1;
    int wsize_ = 10;
};

#endif
