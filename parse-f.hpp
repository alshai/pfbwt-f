#ifndef PARSE_HPP
#define PARSE_HPP

#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cassert>
#include <zlib.h>
#include <sdsl/int_vector.hpp>
#include "hash.hpp"
extern "C" {
#include "utils.h"
#include "gsa/gsacak.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
}

namespace pfbwtf {

struct ParserParams {
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
    Freq(T x) : n(x) {}
    T n = 0;
    T r = 0;
};

struct ntab_entry {
    size_t pos = 0;
    size_t l = 0;
    void clear() {
        pos = 0; l = 0;
    }
};

template <typename Hasher=WangHash>
struct Parser {

    public:

    using UIntType = uint_t;
    using FreqMap = std::map<std::string, Freq<UIntType>>;

    Parser(ParserParams p) : params_(p) {
        check_w(p.w);
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
        uint64_t nseqs(0);
#if !M64
        uint64_t total_l(0);
#endif
        char c('A'), pc('A');
        ntab_entry ne;
        std::string phrase;
        phrase.append(1, Dollar);
        Hasher hf(params_.w);
        while (( l = kseq_read(seq) ) >= 0) {
            if (params_.store_docs) {
                doc_starts_.push_back(pos_);
                doc_names_.push_back(seq->name.s);
            }
#if !M64
            if (total_l + l > 0xFFFFFFFF) {
                fprintf(stderr, "size: %lu\n", total_l + l);
                die("input too long, please use 64-bit version");
            }
            total_l += l;
#endif
            for (size_t i = 0; i < seq->seq.l; ++i) {
                c = std::toupper(seq->seq.s[i]);
                if (params_.trim_non_acgt) {
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
                        continue; // make sure that rest of loop is skipped
                    }
                    if (y < 4 && x > 3) { // record if [^ACGT] run ended
                        ntab_.push_back(ne);
                        ne.clear();
                    }
                } else if (params_.non_acgt_to_a && seq_nt4_table[static_cast<size_t>(c)] > 3) {
                    c = 'A';
                }
                phrase.append(1, c);
                hf.update(c);
                if (hf.hashvalue() % params_.p == 0) {
                    process_phrase(phrase);
                    if (params_.get_sai) {
                        sai.push_back(pos_+1);
                    }
                    phrase.erase(0, phrase.size()-params_.w);
                }
                ++pos_;
                pc = c;
            }
            // record last [^ACGT] run if applicable
            if (params_.trim_non_acgt && seq_nt4_table[static_cast<size_t>(c)] > 3) {
                ntab_.push_back(ne);
            }
            ++nseqs;
        }
        phrase.append(params_.w, Dollar);
        process_phrase(phrase);
        if (params_.get_sai) {
            sai.push_back(pos_ + params_.w);
        }
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

    const std::vector<const char*>& get_parse() const {
        return parse_;
    }

    const std::vector<int_text>& get_parse_ranks() const {
        if (!parse_ranks_.size()) {
            die("parse ranks have not been generated!");
        }
        return parse_ranks_;
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
        if (params_.get_sai) bwsai.push_back(sai[n-1]);
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
                if (params_.get_sai) bwsai.push_back(sai[SA[i]-1]);
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

    const std::vector<ntab_entry>&  get_ntab() const {
        return ntab_;
    }

    const std::vector<UIntType> get_occs() {
        std::vector<UIntType> occs;
        occs.reserve(sorted_phrases_.size());
        for (auto phrase: sorted_phrases_) {
            auto wf = freqs.find(phrase);
            if (wf == freqs.end()) die("there's a problem");
            occs.push_back(wf->second.n);
        }
        return occs;
    }

    const std::vector<const char*>& get_sorted_dict() {
        return sorted_phrases_;
    }

    // sort dictionary, update ranks
    // call when done processing all files
    void finalize() {
        // TODO: add final w Dollars to the parse
        sorted_phrases_.clear();
        sorted_phrases_.reserve(freqs.size());
        for (auto it = freqs.begin(); it != freqs.end(); ++it) {
            sorted_phrases_.push_back(it->first.data());
        }
        std::sort(sorted_phrases_.begin(), sorted_phrases_.end(),
                [](const char* l, const char* r) { return strcmp(l, r) <= 0; });
        size_t rank = 1;
        for (auto phrase: sorted_phrases_) {
            auto& wf = freqs.at(phrase);
            wf.r = rank++;
        }
        parse_ranks_.clear();
        parse_ranks_.reserve(parse_.size());
        for (auto phrase: parse_) {
            auto wf = freqs.at(std::string(phrase));
            parse_ranks_.push_back(wf.r);
        }
    }

    private:

    void inline process_phrase(const std::string& phrase) {
        auto ret = freqs.insert({phrase, Freq<UIntType>(1)});
        if (!ret.second) ret.first->second.n += 1;
        parse_.push_back(ret.first->first.data());
        last_.push_back(phrase[phrase.size()-params_.w-1]);
    }


    FreqMap freqs;
    std::vector<const char*> parse_;
    std::vector<int_text> parse_ranks_;
    std::vector<const char*> sorted_phrases_;
    std::vector<char> last_;
    std::vector<UIntType> sai;
    std::vector<UIntType> doc_starts_;
    std::vector<std::string> doc_names_;
    std::vector<ntab_entry> ntab_;
    ParserParams params_;
    UIntType pos_ = 0;
};
}; // namespace end

#endif
