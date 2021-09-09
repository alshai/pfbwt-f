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
    Marker(size_t p, size_t r, size_t a, size_t s) : textpos(p), refpos(r), allele(a), seqid(s) {}
    size_t textpos;
    size_t refpos;
    size_t allele;
    size_t seqid;
};

class MarkerPositionsWriter {

    public:

    MarkerPositionsWriter() {
    }

    MarkerPositionsWriter(size_t w, FILE* ofp, FILE* log=NULL)
    : wsize_(w)
    , fp_(ofp)
    , log_(log) {}

    ~MarkerPositionsWriter() {
    }

    void update(size_t pos, size_t recpos, int gt, int seqid) {
        if (seqid == -1) {
            fprintf(stderr, "%s ERROR: seqid==-1 not allowed!\n", __FUNCTION__);
        }
        if (seqid_ != -1 && seqid_ != seqid) {
            fprintf(stderr, "%s sequence changed without calling finish_sequence()!\n", __FUNCTION__);
            exit(1);
        }
        // check if new marker is outside window
        while (marker_queue_.size() && marker_queue_.front().textpos + wsize_ <= pos) {
            process_run();
            marker_queue_.pop_front();
        }
        marker_queue_.push_back(Marker(pos, recpos, gt, seqid));
        assert(marker_queue_.back().textpos < marker_queue_.front().textpos + wsize_);
        seqid_ = seqid;
    }

    void finish_sequence() {
        process_run();
        marker_queue_.clear();
        if (markers_to_write_.size()) {
            fwrite(&range_[0], sizeof(uint64_t), 1, fp_);
            fwrite(&range_[1], sizeof(uint64_t), 1, fp_);
            fwrite(markers_to_write_.data(), sizeof(MarkerT), markers_to_write_.size(), fp_);
            fwrite(&delim_, sizeof(uint64_t), 1, fp_);
        }
        markers_to_write_.clear();
        range_[0] = 0; range_[1] = 0;
        seqid_ = -1;
    }

    private:

    void process_run() {
        /* invariant: marker_queue.end() < marker_queue_.begin() + wsize_ */
        uint64_t end;
        // assert(tpos_ <= marker_queue_.front().textpos);  // I don't think this should ever happen
        if (tpos_ + wsize_ <= marker_queue_.front().textpos) { // skip ahead if necessary
            tpos_ = marker_queue_.front().textpos - wsize_ + 1;
        }
        for (auto it = marker_queue_.begin(); it != marker_queue_.end(); ++it) {
            if (!(tpos_ + wsize_ > it->textpos)) {
                end = it->textpos - wsize_;
                assert (end >= tpos_);
                write_markers(tpos_, end, it); // write up to this iter
                tpos_ = end + 1;
            }
        }
        end = marker_queue_.front().textpos;
        write_markers(tpos_, end, marker_queue_.end());
        tpos_ = end + 1;
    }

    void write_markers(uint64_t start, uint64_t end, std::deque<Marker>::iterator it_end) {
        std::vector<MarkerT> markers;
        MarkerT px = -1;
        for (auto it=marker_queue_.begin(); it != it_end; ++it) {
            MarkerT x = create_marker_t(it->refpos, it->allele, it->seqid);
            if (x != px) markers.push_back(x);
            px = x;
        }
        if (start == range_[1] + 1 && markers == markers_to_write_) { // don't write yet if both the same
            range_[1] = end;
        } else {
            if (markers_to_write_.size()) {
                fwrite(&range_[0], sizeof(uint64_t), 1, fp_);
                fwrite(&range_[1], sizeof(uint64_t), 1, fp_);
                fwrite(markers_to_write_.data(), sizeof(MarkerT), markers_to_write_.size(), fp_);
                // fwrite(&x, sizeof(MarkerT), 1, fp_);
                fwrite(&delim_, sizeof(uint64_t), 1, fp_);
            }
            range_[0] = start; range_[1] = end;
            markers_to_write_ = markers;
        }
    }

    uint64_t delim_ = -1;
    size_t wsize_ = 1;
    int seqid_ = -1;
    size_t tpos_ = 0;
    std::deque<Marker> marker_queue_;
    std::vector<MarkerT> markers_to_write_;
    uint64_t range_[2] = {0,0};
    FILE* fp_;
    FILE* log_;
};

template<template<typename> typename ReadConType=VecFileSource>
using MarkerPositions = rle_window_arr<ReadConType>;

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
