#ifndef MARKER_INDEX_HPP
#define MARKER_INDEX_HPP

#include <deque>
#include <cstdio>
#include <cinttypes>
#include "parallel_hashmap/phmap.h"
#include "file_wrappers.hpp"

// false if not equal, true if equal
template<typename T>
bool vec_eq(std::vector<T> a, std::vector<T> b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

struct Marker {
    Marker() {}
    Marker(size_t p, size_t r, size_t a) : textpos(p), refpos(r), allele(a) {}
    size_t textpos;
    size_t refpos;
    size_t allele;
};

class MarkerWindow : public std::deque<Marker> {

    using MDeq = std::deque<Marker>;
    public:

    using iterator = MDeq::iterator;

    MarkerWindow() {}
    MarkerWindow(size_t w) : w_(w) {}

    // returns end() if last element is still within window
    // otherwise, return iterator to last element within begin()->textpos+w
    iterator end_of_window() {
        size_t spos = MDeq::begin()->textpos;
        for (auto it = MDeq::begin(); it != MDeq::end() ; ++it) {
            if (it->textpos >= spos + w_) {
                return it;
            }
        }
        return MDeq::end();
    }

    void push_back(size_t p, size_t r, size_t a) {
        std::deque<Marker>::push_back(Marker(p,r,a));
    }

    size_t get_pos() const { return pos_; }
    void set_pos(size_t p) { pos_ = p; }
    int get_w() const { return w_; }
    void set_w(int w) { w_ = w; }

    private:

    size_t pos_ = 0;
    int w_ = 10;
};


class MarkerIndexWriter {

    public:

    MarkerIndexWriter(int w, FILE* ofp, FILE* log)
        : win_(w)
        , w_(w) { 
            ofp_ = ofp;
            log_ = log;
        }

    ~MarkerIndexWriter() { fprintf(stderr, "n positions w/ windows: %lu\n", nwindows_); }

    void update(size_t rpos, int32_t gt, size_t hpos) {
        if (gt > -1) win_.push_back(hpos, rpos, gt);
        size_t pos = win_.get_pos();
        // skip pos ahead if needed
        if (pos < win_.front().textpos - win_.get_w() + 1) { 
            process_run(); // automatically instigates end of run
            pos = win_.front().textpos - win_.get_w() + 1;
            win_.set_pos(pos);
        }
        // figure out if window is completed
        if (gt > -1) {
            typename MarkerWindow::iterator win_end;
            while ((win_end = win_.end_of_window()) != win_.end()) {
                // process window every position upto and including pos
                for (; pos <= win_.front().textpos; ++pos) {
                    marker_vals.clear();
                    size_t nmarkers = 0;
                    win_.set_pos(pos);
                    for (auto it = win_.begin(); it != win_end; ++it) {
                        if (pos <= it->textpos && it->textpos < pos + win_.get_w()) {
                            ++nmarkers;
                            marker_vals.push_back(it->refpos);
                        } else break;
                    }
                    if (nmarkers) ++nwindows_;
                    if (!vec_eq(pmarker_vals, marker_vals)) {
                        process_run();
                    } 
                    marker_locs.push_back(pos);
                    pmarker_vals = marker_vals;
                }
                win_.pop_front();
            }
        } else { // last rec in file
            for (; pos <= win_.front().textpos; ++pos) {
                marker_vals.clear();
                int nmarkers = 0;
                win_.set_pos(pos);
                for (auto it = win_.begin(); it != win_.end(); ++it) {
                    if (pos <= it->textpos && it->textpos < pos + win_.get_w()) {
                        ++nmarkers;
                        marker_vals.push_back(it->refpos);
                    } else break;
                }
                if (!vec_eq(pmarker_vals, marker_vals)) {
                    process_run();
                    pmarker_vals = marker_vals;
                } 
                marker_locs.push_back(pos);
                process_run();
                if (nmarkers) ++nwindows_;
            }
        }
        ppos = hpos;
    }

    private:

    void process_run() {
        if (marker_locs.size() && pmarker_vals.size()) {
            for (auto l: marker_locs) {
                fwrite(&l, sizeof(uint64_t), 1, ofp_);
            } 
            fwrite(&delim_, sizeof(uint64_t), 1, ofp_);
            for (auto v: pmarker_vals) {
                fwrite(&v, sizeof(uint64_t), 1, ofp_);
            } fwrite(&delim_, sizeof(uint64_t), 1, ofp_);
        }
        marker_locs.clear();
    }

    MarkerWindow win_;
    FILE* ofp_, *log_;
    int w_;

    int32_t ppos;
    uint64_t delim_ = -1;
    size_t nwindows_ = 0;

    std::vector<uint64_t> marker_locs;
    std::vector<uint64_t> marker_vals, pmarker_vals; // switch from uint64_t to std::pair
};

// TODO: many keys have shared values. figure out how to represent this better in a hash table
template<template<typename...>  typename HashMap=phmap::flat_hash_map,
         template<typename>     typename ReadConType=VecFileSource>
class MarkerIndex : public HashMap<uint64_t, std::vector<uint64_t>> {
    
    public:

    MarkerIndex() { }

    MarkerIndex(std::string fname) {
        ReadConType<uint64_t> in_arr(fname);
        std::vector<uint64_t> keys, values;
        int state = 0; // 0 = keys, 1 = values
        for (size_t i = 0; i < in_arr.size(); ++i) {
            if (i == delim_) {
                if (state) { // clear keys, add values to each thing
                    for (auto k: keys) {
                        (*this)[k] = values;
                    }
                    keys.clear();
                    values.clear();
                }
                state = !state;
            }
            else if (state) { // add to values
                values.push_back(in_arr[i]);
            } else { // add to keys
                keys.push_back(in_arr[i]);
            }
        }
    }

    private:

    uint64_t delim_ = -1;

};

#endif
