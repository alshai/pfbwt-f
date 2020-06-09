#ifndef MARKER_INDEX_HPP
#define MARKER_INDEX_HPP

#include <deque>
#include <cstdio>
#include <cinttypes>
#include "file_wrappers.hpp"
#include "rle_window_array.hpp"
#include "marker.hpp"

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


class MarkerPositionsWriter {

    public:

    MarkerPositionsWriter() {}

    MarkerPositionsWriter(int w, FILE* ofp, FILE* log=NULL)
        : win_(w)
        , w_(w) {
            ofp_ = ofp;
            log_ = log;
        }

    ~MarkerPositionsWriter() { }

    void update(size_t rpos, int32_t gt, size_t hpos) {
        if (gt > -1) win_.push_back(hpos, rpos, gt);
        size_t pos = win_.get_pos();
        // skip pos ahead if needed
        if (pos < win_.front().textpos - w_ + 1) {
            process_run(); // automatically instigates end of run
            pos = win_.front().textpos - w_ + 1;
            win_.set_pos(pos);
        }
        if (gt > -1) {
            for (; pos + w_ > win_.front().textpos; ++pos) {
                marker_vals.clear();
                size_t nmarkers = 0;
                typename MarkerWindow::iterator it;
                for (it = win_.begin(); it != win_.end(); ++it) {
                    if (pos + w_ > it->textpos) { // always pos <= m.textpos
                        marker_vals.push_back(create_marker_t(it->refpos, it->allele)); // pack pos and allele
                        ++nmarkers;
                    } else break;
                }
                if (it == win_.end()) { // ignore. more markers to come, probably, with next Marker
                    nmarkers = 0;
                    win_.set_pos(pos);
                    break; // break here, then?
                }
                if (nmarkers) ++nwindows_;
                if (!vec_eq(pmarker_vals, marker_vals)) {
                    process_run();
                }
                marker_locs.push_back(pos);
                pmarker_vals = marker_vals;
                if (pos+1 > win_.front().textpos) {
                    win_.pop_front(); // prune window if necessary
                }
            }
        } else { // last rec in file
            for (; pos <= win_.front().textpos; ++pos) {
                marker_vals.clear();
                int nmarkers = 0;
                win_.set_pos(pos);
                for (auto it = win_.begin(); it != win_.end(); ++it) {
                    if (pos <= it->textpos && it->textpos < pos + win_.get_w()) {
                        ++nmarkers;
                        marker_vals.push_back(create_marker_t(it->refpos, it->allele)); // pack pos and allele
                    } else break;
                }
                if (!vec_eq(pmarker_vals, marker_vals)) {
                    process_run();
                    pmarker_vals = marker_vals;
                }
                marker_locs.push_back(pos);
                // process_run();
                if (nmarkers) ++nwindows_;
            }
            process_run();
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
                fwrite(&v, sizeof(MarkerT), 1, ofp_);
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
    std::vector<MarkerT> marker_vals, pmarker_vals;
};


template<template<typename> typename ReadConType=VecFileSource>
using MarkerPositions = rle_window_arr<ReadConType>;
/*
template<template<typename> typename ReadConType=VecFileSource>
class MarkerPositions : public rle_window_arr<ReadConType> {

    public:

    MarkerPositions() {}
    MarkerPositions(std::string fname) : rle_window_arr<ReadConType>(fname) {}

    bool has_markers(uint64_t i) const {
        return this->has_entry(i);
    }

    std::vector<uint64_t> get_markers(uint64_t i) const {
        return this->at(i);
    }
};
*/

template <typename MPos=MarkerPositions<>>
void write_marker_array(std::string mai_fname, std::string sa_fname, std::string output = "") {
    FILE* sa_fp = sa_fname == "-" ? stdin : fopen(sa_fname.data(), "rb");
    FILE* ofp = fopen(output == "" ? "out" : output.data(), "wb");
    MPos mai(mai_fname);
    constexpr uint64_t delim = -1;
    uint64_t s;
    uint64_t i = 0;
    std::vector<uint64_t> markers, pmarkers, locs;
    while (fread(&s, sizeof(uint64_t), 1, sa_fp) == 1) {
        markers.clear();
        markers = mai.at(s);
        if (!vec_eq(markers, pmarkers)) {
            if (pmarkers.size()) {
                fwrite(&locs.front(), sizeof(uint64_t), 1, ofp);
                fwrite(&locs.back(), sizeof(uint64_t), 1, ofp);
                for (auto m: pmarkers) {
                    fwrite(&m, sizeof(uint64_t), 1, ofp);
                }
                fwrite(&delim, sizeof(uint64_t), 1, ofp);
            }
            locs.clear();
        }
        locs.push_back(i);
        pmarkers = markers;
        ++i;
    }
    if (pmarkers.size()) {
        fwrite(&locs.front(), sizeof(uint64_t), 1, ofp);
        fwrite(&locs.back(), sizeof(uint64_t), 1, ofp);
        for (auto m: pmarkers) {
            fwrite(&m, sizeof(uint64_t), 1, ofp);
        }
        fwrite(&delim, sizeof(uint64_t), 1, ofp);
    }
    fclose(sa_fp);
    fclose(ofp);
}

template<template<typename> typename ReadConType=VecFileSource>
using MarkerArray = rle_window_arr<ReadConType>;

#endif
