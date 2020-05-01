#ifndef MARKER_INDEX_HPP
#define MARKER_INDEX_HPP

#include <deque>
#include <cstdio>

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

    void update(size_t rpos, int32_t gt, size_t hpos) {
        if (gt > -1) win_.push_back(hpos, rpos, gt);
        size_t pos = win_.get_pos();
        // skip pos ahead if needed
        if (pos < win_.front().textpos - win_.get_w() + 1) {
            pos = win_.front().textpos - win_.get_w() + 1;
            win_.set_pos(pos);
        }
        // figure out if window is completed
        if (gt > -1) {
            typename MarkerWindow::iterator win_end;
            while ((win_end = win_.end_of_window()) != win_.end()) {
                // process window every position upto and including pos
                for (; pos <= win_.front().textpos; ++pos) {
                    win_.set_pos(pos);
                    fwrite(&pos, sizeof(size_t), 1, ofp_);
                    fprintf(log_, "%lu: ", pos);
                    for (auto it = win_.begin(); it != win_end; ++it) {
                        if (pos <= it->textpos && it->textpos < pos + win_.get_w()) {
                            fprintf(log_, "(%lu %lu %lu) ", it->textpos, it->refpos, it->allele);
                            fwrite(&(it->refpos), sizeof(size_t), 1, ofp_);
                            fwrite(&(it->allele), sizeof(size_t), 1, ofp_);
                        } else break;
                    }
                    fwrite(&negone, sizeof(size_t), 1, ofp_);
                    fprintf(log_, "\n");
                }
                win_.pop_front();
            }
        } else { // last rec in file
            for (; pos <= win_.front().textpos; ++pos) {
                win_.set_pos(pos);
                fwrite(&pos, sizeof(size_t), 1, ofp_);
                fprintf(log_, "%lu: ", pos);
                for (auto it = win_.begin(); it != win_.end(); ++it) {
                    if (pos <= it->textpos && it->textpos < pos + win_.get_w()) {
                        fprintf(log_, "(%lu %lu %lu) ", it->textpos, it->refpos, it->allele);
                        fwrite(&(it->refpos), sizeof(size_t), 1, ofp_);
                        fwrite(&(it->allele), sizeof(size_t), 1, ofp_);
                    } else break;
                }
                fwrite(&negone, sizeof(size_t), 1, ofp_);
                fprintf(log_, "\n");
            }
        }
    }

    private:

    MarkerWindow win_;
    FILE* ofp_, *log_;
    int w_;

    size_t zero = 0;
    size_t negone = -1;
};

#endif
