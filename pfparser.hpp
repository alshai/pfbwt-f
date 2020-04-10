#ifndef PFPARSE_HPP
#define PFPARSE_HPP

/* Author: Taher Mun
 * Date  : April 10, 2020
 */

#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <cassert>
#include <zlib.h>
#include "hash.hpp"
extern "C" {
#include "utils.h"
#include "gsa/gsacak.h"
#ifndef AC_KSEQ_H
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
}

namespace pfbwtf {

struct PfParserParams {
    size_t w = 10;
    size_t p = 100;
    bool get_sai = false;
    bool store_docs = false;
    bool verbose = false;
    // if bothm trim_non_acgt and non_acgt_to_a are set, trim_non_acgt takes precedence
    bool trim_non_acgt = false;
    bool non_acgt_to_a = false;
};

template<typename T>
struct Freq {
    Freq() {}
    Freq(T x) : n(x) {}
    Freq(T x, T y) : n(x), r(y) {}
    T n = 0;
    T r = 0;
    bool operator==(const Freq<T>& rhs) const {
        return (rhs.n == n && rhs.r == r);
    }
    bool operator!=(const Freq<T>& rhs) const {
        return (rhs.n != n || rhs.r != r);
    }
};

struct ntab_entry {
    size_t pos = 0;
    size_t l = 0;
    void clear() {
        pos = 0; l = 0;
    }
};

template<typename UIntType>
using FreqMap = std::map<std::string, Freq<UIntType>>;

template <typename Hasher=WangHash>
struct PfParser {

    public:

    using UIntType = uint_t;
    using IntType = int_text;

    PfParser() {}

    PfParser(PfParserParams p) : params_(p) {
        check_w(p.w);
    }

    // load from serialized dict, occs and parse ranks, 
    // and others as user sees fit
    // assuming dict is sorted, occs same order as dict
    PfParser(PfParserParams params, 
           std::vector<std::string>& sorted_phrases, 
           std::vector<IntType>&& parse_ranks,
           std::vector<UIntType>&& doc_starts = std::vector<UIntType>(),
           std::vector<std::string>&& doc_names = std::vector<std::string>(),
           std::vector<ntab_entry>&& ntab = std::vector<ntab_entry>()
           ) 
        : parse_ranks_(parse_ranks)
        , doc_starts_(doc_starts)
        , doc_names_(doc_names)
        , ntab_(ntab)
        , params_(params)
    {
        sorted_phrases_.reserve(sorted_phrases.size());
        parse_.reserve(parse_ranks_.size());
        last_.reserve(parse_ranks_.size());
        std::string phrase;
        for (size_t i = 0; i < parse_ranks_.size()-1; ++i) {
            phrase.assign(sorted_phrases[parse_ranks_[i]-1]);
            pos_ += pos_ ? phrase.size() - params_.w : phrase.size() - 1;
            process_phrase(phrase);
        }
        last_phrase_ = sorted_phrases[parse_ranks_.back()-1];
        // account for first w characters AND last w Dollars here (hence "2*")
        pos_ += last_phrase_.size() - (2 * params_.w) + 1; 
        last_phrase_.erase(last_phrase_.size() - params_.w, params_.w);
        finalize();
    }

    PfParser(const PfParser& rhs) 
        : freqs_(rhs.freqs_)
        , parse_ranks_(rhs.parse_ranks_)
        , last_(rhs.last_)
        , sai_(rhs.sai_)
        , doc_starts_(rhs.doc_starts_)
        , doc_names_(rhs.doc_names_)
        , params_(rhs.params_)
        , pos_(rhs.pos_)
        , nseqs_(rhs.nseqs_)
    {
        sort_dict();
        regenerate_parse();
        last_phrase_ = std::string(parse_.back());
    }

    PfParser(PfParser&& rhs) 
        : freqs_(std::move(rhs.freqs_))
        , parse_ranks_(std::move(rhs.parse_ranks_))
        , last_(std::move(rhs.last_))
        , sai_(std::move(rhs.sai_))
        , doc_starts_(std::move(rhs.doc_starts_))
        , doc_names_(std::move(rhs.doc_names_))
        , params_(std::move(rhs.params_))
        , pos_(std::move(rhs.pos_))
        , nseqs_(std::move(rhs.nseqs_))
    {
        sort_dict();
        regenerate_parse();
        last_phrase_ = std::string(parse_.back());
    }

    PfParser& operator=(const PfParser& rhs) {
        params_ = rhs.params_;
        freqs_ = rhs.freqs_;
        parse_ranks_ = rhs.parse_ranks_;
        last_ = rhs.last_;
        sai_ = rhs.sai_;
        pos_ = rhs.pos_;
        doc_starts_ = rhs.doc_starts_;
        doc_names_ = rhs.doc_names_;
        nseqs_ = rhs.nseqs_;
        // parse_ and sorted_phrases_ are pointers to freq_ keys, so they need to be regenerated just in case
        sort_dict();
        regenerate_parse();
        last_phrase_ = std::string(parse_.back());
        return *this;
    }

    PfParser& operator=(PfParser&& rhs) {
        params_ = std::move(rhs.params_);
        freqs_ = std::move(rhs.freqs_);
        parse_ranks_ = std::move(rhs.parse_ranks_);
        last_ = std::move(rhs.last_);
        sai_ = std::move(rhs.sai_);
        pos_ = rhs.pos_;
        doc_starts_ = std::move(rhs.doc_starts_);
        doc_names_ = std::move(rhs.doc_names_);
        nseqs_ = rhs.nseqs_;
        // parse_ and sorted_phrases_ are pointers to freq_ keys, so they need to be regenerated just in case
        sort_dict();
        regenerate_parse();
        return *this;
    }

    // append information from another parse 
    // remember to use .finalize() after finishing using +=!
    PfParser& operator+=(const PfParser& rhs) {
        if (!freqs_.size()) return operator=(rhs);
        if (rhs.params_.w != params_.w) exit(1);
        if (rhs.params_.p != params_.p) exit(1);
        // retrieve final phrase of this parse 
        std::string phrase(parse_.back());
        // decrement count of last phrase in frequency table
        auto it = freqs_.find(phrase);
        if (it == freqs_.end()) die("last phrase was supposed to be in dict");
        if (it->second.n) {
            --(it->second.n);
        }
        // made this an 'if' instead of an 'else' on purpose
        if (!it->second.n) {
            freqs_.erase(it);
        }
        // erase dollars if present
        if (phrase.back() == Dollar) {
            phrase.erase(phrase.size() - params_.w, params_.w);
            pos_ -= params_.w;
        }
        // update other data structures as well to reflect removal of last phease
        parse_.pop_back();
        last_.pop_back();
        if (params_.get_sai) sai_.pop_back();
        // parse_ranks_.pop_back(); // don't really need to do this here bc ranks will be regenerated later
        // TODO: deal with doc_starts_, doc_names_, ntab_
        if (rhs.parse_[0][0] != Dollar) die("rhs parser malformed");
        // join last phrase of this parse and first phrase of next parse
        // first, load hasher with `w` As
        Hasher hf(params_.w);
        for (size_t i = 0; i < params_.w; ++i) hf.update('A');
        char c; // , pc;
        // we want to look at the last four As and the first four of the window
        // (ie. from AAAA_ to A____)
        // because we're removing the dollar, we also have to redo _____
        std::string window(rhs.parse_[0]+1, params_.w);
        assert(window[0] != Dollar);
        for (size_t i = 0; i < window.size(); ++i) {
            c = window[i];
            phrase.append(1, c);
            hf.update(c);
            if (hf.hashvalue() % params_.p == 0) {
                ++pos_; // TODO: make this less hacky
                process_phrase(phrase);
                phrase.erase(0, phrase.size()-params_.w);
                --pos_; // TODO: here too
            }
            ++pos_;
        }
        // at this point last_phrase ends with first four characters of rhs
        // we know that rhs.parse_[0] already ends in a window so just join it to last_phrase
        phrase.append(rhs.parse_[0]+params_.w+1);
        pos_ += strlen(rhs.parse_[0])-params_.w-1;
        process_phrase(phrase);
        // concatenate the rest of rhs.parse_
        parse_.reserve(parse_.size() + rhs.parse_.size()-1);
        for (size_t i = 1; i < rhs.parse_.size(); ++i) {
            phrase.assign(rhs.parse_[i]);
            pos_ += phrase.size() - params_.w;
            process_phrase(phrase);
        } // NOTE: if rhs is finalized, there are also Dollars at end
        last_phrase_ = parse_.back();
        // TODO: concatenate doc_names
        //       this means add pos_ to each doc
        nseqs_ += rhs.nseqs_;
        return *this;
    }

    bool operator==(const PfParser& rhs) {
        if (pos_ != rhs.pos_) return false;
        if (parse_ranks_.size() != rhs.parse_ranks_.size()) return false;
        if (last_.size() != rhs.last_.size()) return false;
        if (freqs_.size() != rhs.freqs_.size()) return false;
        if (params_.get_sai) {
            if (sai_.size() != rhs.sai_.size()) return false;
        }
        for (const auto& pair: freqs_) {
            const auto it2 = rhs.freqs_.find(pair.first);
            if (it2 == rhs.freqs_.end()) return false;
            if (it2->second != pair.second) return false;
        }
        for (size_t i = 0; i < parse_ranks_.size(); ++i) {
            if (parse_ranks_[i] != rhs.parse_ranks_[i]) return false;
        }
        for (size_t i = 0; i < last_.size(); ++i) {
            if (last_[i] != rhs.last_[i]) return false;
        }
        for (size_t i = 0; i < parse_.size(); ++i) {
            if (strcmp(parse_[i], rhs.parse_[i])) return false;
        }
        for (size_t i = 0; i < sorted_phrases_.size(); ++i) {
            if (strcmp(sorted_phrases_[i], rhs.sorted_phrases_[i])) return false;
        }
        if (params_.get_sai) {
            for (size_t i = 0; i < sai_.size(); ++i) {
                if (sai_[i] != rhs.sai_[i]) return false;
            }
        }
        return true;
    }

    // stores parse information from a fasta file
    size_t add_fasta(std::string fasta_fname) {
        gzFile fp;
        if (fasta_fname == "-") {
            fp = gzdopen(fileno(stdin), "r");
        } else {
            fp = gzopen(fasta_fname.data(), "r");
        }
        if (fp == NULL) die("failed to open file!\n");
        kseq_t* seq = kseq_init(fp);
        int l;
#if !M64
        uint64_t total_l(0);
#endif
        ntab_entry ne;
        char c('A'), pc('A');
        std::string phrase(last_phrase_);
        if (!pos_) {
            phrase.append(1, Dollar);
            ++pos_;
        }
        Hasher hf(params_.w);
        while (( l = kseq_read(seq) ) >= 0) {
            nseqs_ += 1;
            if (params_.store_docs) {
                // subtract 1 Dollar and w*(npreviousdocs) Dollars
                doc_starts_.push_back(pos_ - 1 - (params_.w * doc_starts_.size()));
                doc_names_.push_back(seq->name.s);
            }
#if !M64
            if (total_l + l > 0xFFFFFFFF) {
                fprintf(stderr, "size: %lu\n", total_l + l);
                die("input too long, please use 64-bit version");
            }
            total_l += l;
#endif
            // loop through incoming sequence, and add `w` extra 'A's to
            // sequence so we can add additional sequences later
            for (size_t i = 0; i < seq->seq.l + params_.w; ++i) {
                // this if-else ensures `w` As will be appended to end of seq
                c = i < seq->seq.l ? std::toupper(seq->seq.s[i]) : 'A';
                if (params_.trim_non_acgt) {
                    if (update_ntab(pc, c, ne)) continue;
                } else if (params_.non_acgt_to_a && seq_nt4_table[static_cast<size_t>(c)] > 3) {
                    c = 'A';
                }
                phrase.append(1, c);
                hf.update(c);
                if (pos_ > params_.w && hf.hashvalue() % params_.p == 0) {
                    process_phrase(phrase);
                    phrase.erase(0, phrase.size()-params_.w);
                }
                ++pos_;
                pc = c;
            }
            // record last [^ACGT] run if applicable
            if (params_.trim_non_acgt && seq_nt4_table[static_cast<size_t>(c)] > 3) {
                ntab_.push_back(ne);
            }
        }
        last_phrase_ = phrase;
        // phrase.append(params_.w, Dollar);
        // process_phrase(phrase);
        //
        kseq_destroy(seq);
        gzclose(fp);
        return pos_;
    }

    void check_w(size_t x) {
        if (x > 32){
            fprintf(stderr, "window size w must be < 32!\n");
            exit(1);
        }
    }

    // generates bwlast and ilist (and bwsai)
    template<typename OutFn>
    void bwt_of_parse(OutFn out_fn) {
        auto occs = get_occs();
        // these will get passed to out_fn at end
        std::vector<char> bwlast;
        std::vector<UIntType> ilist;
        std::vector<UIntType> bwsai;

        size_t n; // size of parse_ranks, minus the last EOS character
        // TODO: support large parse sizes
        if (!parse_ranks_.size()) get_parse_ranks();
        if (parse_ranks_.size() == 1) {
            die("error: only one dict word total. Re-run with a smaller p modulus");
        }
#if !M64
        // if in 32 bit mode, the number of words is at most 2^31-2
        if(parse_ranks_.size() > 0x7FFFFFFE) {
            fprintf(stderr, "parse ranks size: %lu\n", parse_ranks_.size());
            die("Input containing more than 2^31-2 phrases! Please use 64 bit version");
        }
#else
        // if in 64 bit mode, the number of words is at most 2^32-2 (for now)
        if(parse_ranks_.size() > 0xFFFFFFFEu) {
            fprintf(stderr, "parse ranks size: %lu\n", parse_ranks_.size());
            die("Input containing more than 2^32-2 phrases! This is currently a hard limit");
        }
#endif
        // add EOS if it's not already there
        if (parse_ranks_[parse_ranks_.size()-1]) {
            n = parse_ranks_.size();
            parse_ranks_.push_back(0);
        } else n = parse_ranks_.size() - 1;
        // TODO: calculate k
        size_t k = 0;
        for (size_t i = 0; i < n; ++i) {
            k = parse_ranks_[i] > k ? parse_ranks_[i] : k;
        }
        if (params_.get_sai) {
            bwsai.clear();
            bwsai.reserve(n);
        }
        // compute S.A.
        // we assign instead of reserve in order to be able to use .size()
        assert(sizeof(UIntType) >= sizeof(uint_t));
        std::vector<UIntType> SA(n+1, 0);
        // fprintf(stderr, "Computing S.A. of size %ld over an alphabet of size %ld\n",n+1,k+1);
        int depth = sacak_int(parse_ranks_.data(), SA.data(), n+1, k+1);
        if (depth < 0) die("Error computing SA");
        // if(depth>=0) fprintf(stderr, "S.A. computed with depth: %d\n", depth);
        // else die("Error computing the S.A.");
        // transform S.A. to BWT in place
        assert(SA[0] == n);
        SA[0] = parse_ranks_[n-1];
        bwlast.clear();
        bwlast.reserve(n+1);
        bwlast.push_back(last_[n-2]);
        if (params_.get_sai) bwsai.push_back(sai_[n-1]);
        for (size_t i = 1; i < n+1; ++i) {
            if (!SA[i]) {
                SA[i] = 0;
                bwlast.push_back(0);
                // TODO: write to bwsai
                if (params_.get_sai) bwsai.push_back(0);
            } else {
                if (SA[i] == 1) {
                    bwlast.push_back(last_[n-1]);
                } else {
                    bwlast.push_back(last_[SA[i]-2]);
                }
                if (params_.get_sai) bwsai.push_back(sai_[SA[i]-1]);
                SA[i] = parse_ranks_[SA[i] - 1];
            }
        }
        std::vector<UIntType> F(occs.size()+1, 0);
        F[1] = 1;
        for (size_t i = 2; i < occs.size() + 1; ++i) {
            F[i] = F[i-1] + occs[i-2];
        }
        assert(F[occs.size()] + occs[occs.size()-1] == n+1);
        // TODO: do we want to store ilist as a bitvector directly?
        ilist.resize(n+1, 0);
        for (size_t i = 0; i < n + 1; ++i) {
            ilist[F[SA[i]]++] = i;
        }
        // ilist_processor(ilist);
        assert(ilist[0]==1);
        assert(SA[ilist[0]] == 0);
        out_fn(bwlast, ilist, bwsai);
    }

    size_t get_parse_size() const { return parse_ranks_.size(); }

    const std::vector<UIntType>& get_doc_starts() const {
        return doc_starts_;
    }

    const std::vector<std::string>& get_doc_names() const {
        return doc_names_;
    }

    const std::vector<UIntType> get_occs() {
        std::vector<UIntType> occs;
        occs.reserve(sorted_phrases_.size());
        for (auto phrase: sorted_phrases_) {
            auto wf = freqs_.find(phrase);
            if (wf == freqs_.end()) die("there's a problem");
            occs.push_back(wf->second.n);
        }
        return occs;
    }

    // sort dictionary, update ranks
    // call when done processing all files
    void finalize() {
        if (last_phrase_.back() != Dollar)  {
            last_phrase_.append(params_.w, Dollar);
            pos_ += params_.w - 1; // because we added w Dollars
            process_phrase(last_phrase_);
        }
        sort_dict();
        generate_ranks();
    }

    void sort_dict() {
        sorted_phrases_.clear();
        sorted_phrases_.reserve(freqs_.size());
        for (auto it = freqs_.begin(); it != freqs_.end(); ++it) {
            sorted_phrases_.push_back(it->first.data());
        }
        std::sort(sorted_phrases_.begin(), sorted_phrases_.end(),
                [](const char* l, const char* r) { return strcmp(l, r) <= 0; });
    }

    void generate_ranks() {
        if (!sorted_phrases_.size()) sort_dict();
        size_t rank = 1;
        for (auto phrase: sorted_phrases_) {
            auto& wf = freqs_.at(phrase);
            wf.r = rank++;
        }
        parse_ranks_.clear();
        if (parse_.size()) parse_ranks_.reserve(parse_.size());
        for (auto phrase: parse_) {
            auto wf = freqs_.at(std::string(phrase));
            parse_ranks_.push_back(wf.r);
        }
    }

    void regenerate_parse() {
        if (!sorted_phrases_.size()) { sort_dict(); }
        if (!parse_ranks_.size()) { generate_ranks(); }
        parse_.clear();
        for (auto r: parse_ranks_) {
            parse_.push_back(sorted_phrases_[r - 1]);
        }
    }

    size_t get_n() {  // includes As at end of each seq
        //             Dollars
        return pos_ - params_.w; 
    }

    const std::vector<UIntType>& get_sai() const { return sai_; }
    const std::vector<char>& get_last() const { return last_; }
    const FreqMap<UIntType>& get_freqs() const { return freqs_; }
    const std::vector<const char*>& get_parse() const { return parse_; }
    const std::vector<int_text>& get_parse_ranks() const { return parse_ranks_; }
    size_t get_pos() const { return pos_; }
    const std::vector<const char*>& get_sorted_phrases() const { return sorted_phrases_; }
    const std::vector<ntab_entry>&  get_ntab() const { return ntab_; }
    const PfParserParams get_params() const { return params_; }

    private:

    bool inline update_ntab(char pc, char c, ntab_entry& ne) {
        char x = seq_nt4_table[static_cast<size_t>(pc)];
        char y = seq_nt4_table[static_cast<size_t>(c)];
        if (y > 3) { // skip if nonACGT
            if (x < 4) { // new N run
                ne.l = 1;
                ne.pos = pos_-1;
            } else {
                ne.l += 1;
            }
            pc = c;
            return true;
        }
        if (y < 4 && x > 3) { // record if [^ACGT] run ended
            ntab_.push_back(ne);
            ne.clear();
        }
        return false;
    }

    void inline process_phrase(const std::string& phrase) {
        auto ret = freqs_.insert({phrase, Freq<UIntType>(1)});
        if (!ret.second) ret.first->second.n += 1;
        parse_.push_back(ret.first->first.data());
        last_.push_back(phrase[phrase.size()-params_.w-1]);
        if (params_.get_sai) sai_.push_back(pos_);
    }


    FreqMap<UIntType> freqs_; // # unique words
    std::vector<const char*> parse_; // # words
    std::vector<int_text> parse_ranks_; // # words
    std::vector<const char*> sorted_phrases_; // # unique words
    std::vector<char> last_; // # words
    std::vector<UIntType> sai_; // # words
    std::vector<UIntType> doc_starts_;
    std::vector<std::string> doc_names_;
    std::vector<ntab_entry> ntab_;
    PfParserParams params_;
    UIntType pos_ = 0; // total characters in parse, INCLUDING Dollar chars!
    std::string last_phrase_ = "";
    size_t nseqs_ = 0;
};
}; // namespace end

#endif // PFPARSE_HPP
