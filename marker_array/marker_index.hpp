#ifndef MARKER_INDEX_HPP
#define MARKER_INDEX_HPP

#include <deque>
#include <cstdio>
#include <cinttypes>
#include "parallel_hashmap/phmap.h"
#include "file_wrappers.hpp"
#include "sdsl_bv_wrappers.hpp"

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

    ~MarkerIndexWriter() {
        // fprintf(stderr, "n positions w/ windows: %lu\n", nwindows_);
        // fprintf(stderr, "last position: %d\n", ppos);
    }

    void update(size_t rpos, int32_t gt, size_t hpos) {
        if (gt > -1) win_.push_back(hpos, rpos, gt);
        size_t pos = win_.get_pos();
        // skip pos ahead if needed
        if (pos < win_.front().textpos - w_ + 1) {
            fprintf(log_, "skipping from %lu ", pos);
            process_run(); // automatically instigates end of run
            pos = win_.front().textpos - w_ + 1;
            fprintf(log_, "to %lu \n", pos);
            win_.set_pos(pos);
        } // else if (inc_) ++pos;
        if (gt > -1) {
            for (; pos + w_ > win_.front().textpos; ++pos) {
                fprintf(log_, "%lu: ", pos);
                marker_vals.clear();
                size_t nmarkers = 0;
                typename MarkerWindow::iterator it;
                for (it = win_.begin(); it != win_.end(); ++it) {
                    if (pos + w_ > it->textpos) { // always pos <= m.textpos
                        marker_vals.push_back(it->refpos);
                        ++nmarkers;
                    } else break;
                }
                if (it == win_.end()) { // ignore. more markers to come, probably, with next Marker
                    nmarkers = 0;
                    win_.set_pos(pos);
                    fprintf(log_, " break\n");
                    break; // break here, then?
                }
                if (nmarkers) ++nwindows_;
                if (!vec_eq(pmarker_vals, marker_vals)) {
                    fprintf(log_, "processing... (");
                    for (auto x: marker_locs) {
                        fprintf(log_, "%lu ", x);
                    }
                    fprintf(log_, ") [");
                    for (auto x: pmarker_vals) {
                        fprintf(log_, "%lu ", x);
                    }
                    fprintf(log_, ") ]");
                    process_run();
                    // inc_ = true;
                } // else inc_ = false;
                marker_locs.push_back(pos);
                pmarker_vals = marker_vals;
                if (pos+1 > win_.front().textpos) {
                    fprintf(log_, "(%lu -> %lu) pruning %lu", pos, pos+1, win_.front().refpos);
                    win_.pop_front(); // prune window if necessary
                }
                fprintf(log_, "\n");
            }
            /*
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
                        inc_ = true;
                    } else inc_ = false;
                    marker_locs.push_back(pos);
                    pmarker_vals = marker_vals;
                }
                win_.pop_front();
            }
            */
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
            uint64_t front = marker_locs.front();
            uint64_t back = marker_locs.back();
            fwrite(&front, sizeof(uint64_t), 1, ofp_);
            fwrite(&back, sizeof(uint64_t), 1, ofp_);
            for (auto v: pmarker_vals) {
                fwrite(&v, sizeof(uint64_t), 1, ofp_);
            }
            fwrite(&delim_, sizeof(uint64_t), 1, ofp_);
        }
        marker_locs.clear();
    }

    MarkerWindow win_;
    FILE* ofp_, *log_;
    int w_;

    int32_t ppos;
    uint64_t delim_ = -1;
    size_t nwindows_ = 0;
    bool inc_ = false;

    std::vector<uint64_t> marker_locs;
    std::vector<uint64_t> marker_vals, pmarker_vals, vbuf; // switch from uint64_t to std::pair
};

template<template<typename> typename ReadConType=VecFileSource>
class MarkerIndex {

    public:

    MarkerIndex() { }

    MarkerIndex(std::string fname) {
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
                rs = keys[0];
                re = keys[1];
                marker_windows_.push_back(values);
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
     *    +  -X+  -X+-X+ -XX+ -
     * 000100001000010010000100
     * 000000100001001000100001
     */
    bool has_markers(uint64_t i) {
        return (run_starts_.rank(i+1) == run_ends_.rank(i) + 1);
    }

    std::vector<uint64_t> get_markers(uint64_t i) {
        uint64_t srank = run_starts_.rank(i+1);
        uint64_t erank = run_ends_.rank(i);
        if (srank != erank + 1) return std::vector<uint64_t>();
        else return marker_windows_[srank-1];
    }

    private:

    uint64_t get_last_position(ReadConType<uint64_t>& in_arr) {
        uint64_t i = in_arr.size() - 2; // in_arr[-1] is delim_
        while (in_arr[i] != delim_) {
            --i;
        }
        return in_arr[i+2]; //  delim keys[0] keys[1] val1 val2 ... delim EOF
    }

    uint64_t delim_ = -1;
    bv_rs<> run_starts_;
    bv_rs<> run_ends_;
    std::vector<std::vector<uint64_t>> marker_windows_;
};

// template<template<typename...>  typename HashMap=phmap::flat_hash_map,
//          template<typename>     typename ReadConType=VecFileSource>
// class MarkerIndex : public HashMap<uint64_t, std::vector<uint64_t>> {
// 
//     public:
// 
//     MarkerIndex() { }
// 
//     MarkerIndex(std::string fname) {
//         ReadConType<uint64_t> in_arr(fname);
//         uint64_t keys[2];
//         std::vector<uint64_t> values;
//         int state = 0;
//         for (size_t i = 0; i < in_arr.size(); ++i) {
//             if (in_arr[i] == delim_) {
//                 if (keys[0] == keys[1]) {
//                     (*this)[keys[0]] = values;
//                 } else {
//                     // fprintf(stderr, "processing keys %lu to %lu\n", keys[0], keys[1]);
//                     for (uint64_t k = keys[0]; k <= keys[1]; ++k) {
//                         (*this)[k] = values;
//                     }
//                 }
//                 // keys.clear();
//                 values.clear();
//                 state = 0;
//             }
//             else if (state < 2) {
//                 keys[state++] = in_arr[i];
//             } else {
//                 values.push_back(in_arr[i]);
//             }
//         }
//     }
// 
//     std::vector<uint64_t> get_markers(uint64_t i) {
//         return (*this)[i];
//     }
// 
//     private:
// 
//     uint64_t delim_ = -1;
// 
// };


#endif

