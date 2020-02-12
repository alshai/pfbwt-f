#ifndef PFBWT_HPP
#define PFBWT_HPP

#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>
#include <fcntl.h>
#include "sdsl/bit_vectors.hpp"
#include "sdsl_bv_wrappers.hpp"
extern "C" {
#include <sys/mman.h>
#include "utils.h"
#include "gsa/gsacak.h"
}

namespace pfbwtf {

struct SuffixT {
    SuffixT(uint8_t c, uint_t i, uint_t w) : bwtc(c), bwtp(i), wordi(w) {}
    uint8_t bwtc;
    uint_t bwtp;
    uint_t wordi;
    bool operator<(SuffixT& r);
};

bool SuffixT::operator<(SuffixT& r) {
    return this->bwtp < r.bwtp;
}

enum class RunType {OTHER, START, END};
enum class Difficulty {EASY1, EASY2, HARD};

struct sa_fn_arg {
    sa_fn_arg(uint_t p, uint_t s, RunType r = RunType::OTHER, Difficulty d = Difficulty::EASY1) :
        pos(p), sa(s), run_t(r), dif(d) {}
    uint_t pos;
    uint_t sa;
    RunType run_t;
    Difficulty dif;
};


template<template <typename, typename...> typename ReadConType,
         template <typename, typename...> typename WriteConType
         >
class PrefixFreeBWT {

    public:

    using UIntType = uint_t;
    using IntType = int_t;

    PrefixFreeBWT(std::string prefix, size_t win_size, bool sa = false, bool ssa = false, bool verb = false) :
        fname(prefix),
        w ( win_size),
        dict ( WriteConType<uint8_t>(prefix + "." + EXTDICT)),
        bwlast ( ReadConType<uint8_t>(prefix + "." + EXTBWLST)),
        ilist ( ReadConType<UIntType>(prefix + "." + EXTILIST)),
        build_sa(sa), build_rssa(ssa),
        any_sa(sa | ssa),
        verbose(verb)
    {
        if (verbose) fprintf(stderr, "loaded files\n");
        if (sa && ssa) die("cannot activate both SA and sampled-SA options!");
        dsize = dict.size();
        if (verbose) fprintf(stderr, "creating ilist idx\n");
        if (sa && ssa) die("cannot activate both SA and sampled-SA options!");
        load_ilist_idx(prefix);
        if (sa || ssa) bwsai = ReadConType<UIntType>(prefix + "." + EXTBWSAI);
    }

#define get_word_suflen(i, d, s) \
    d = dict_idx.rank(i); \
    s = d>=dwords ? dsize-i : dict_idx.select(d+1) - i;

#define UPDATE_SA(bwtc, bwtp, d) \
    sa = bwtp - suff_len; \
    if (build_sa) { \
        sa_fn(sa_fn_arg(pos, sa)); \
    } else if (build_rssa && bwtc != pbwtc) { \
        sa_fn(sa_fn_arg(pos,   sa,  RunType::START, d)); \
        sa_fn(sa_fn_arg(pos-1, psa, RunType::END,   d)); \
    } \
    psa = sa;

#define UPDATE_BWT(bwtc) \
    bwt_fn(bwtc); \
    pbwtc = bwtc;

     /* uses LCP of dict to build BWT (less memory, more time)
     */
    template<typename BFn, typename SFn>
    void generate_bwt_lcp(BFn bwt_fn, SFn sa_fn) {
        if (verbose) fprintf(stderr, "generating dict suffixes\n");
        sort_dict_suffixes(true); // build gSA and gLCP of dict
        // start from SA item that's not EndOfWord or EndOfDict
        size_t next, suff_len, wordi;
        uint8_t bwtc, pbwtc = 0;
        size_t easy_cases = 0, hard_cases = 0;
        size_t pos = 0;
        UIntType sa, psa = 1;
        if (verbose) fprintf(stderr, "processing words to build BWT\n");
        for (size_t i = dwords+w+1; i<dsize; i=next) {
            next = i+1;
            get_word_suflen(gsa[i], wordi, suff_len);
            if (suff_len <= w) continue; // ignore small suffixes
            // full word case
            if (gsa[i] == 0 || dict_idx[gsa[i]-1] == 1) {
                auto word_ilist = get_word_ilist(wordi);
                for (auto j: word_ilist) {
                    bwtc = bwlast[j];
                    if (any_sa) {
                        if (wordi > 0) {
                            UPDATE_SA(bwtc, bwsai[j], Difficulty::EASY1);
                        } else if (build_rssa) { // wordi == 0
                            sa = bwsai[0] - this->w;
                            sa_fn(sa_fn_arg(pos, sa, RunType::START, Difficulty::EASY1));
                            psa = sa; // save current sa 
                        }
                        /*
                        if (update_da) {
                            UPDATE_DA(sa);
                        }
                        */
                    }
                    UPDATE_BWT(bwtc);
                    ++pos;
                    ++easy_cases;
                }
            } else { // hard case!
                // look at all the sufs that share LCP[suf]==this_suffixlen
                size_t nwordi, nsuff_len;
                std::vector<uint8_t> chars;
                std::vector<uint64_t> words;
                uint8_t c, pc = gsa[i]-1 ? dict[gsa[i]-1] : 0;
                chars.push_back(pc);
                words.push_back(wordi);
                bool same_char = true;
                size_t j;
                for (j = i + 1; j < dsize && glcp[j] >= (IntType) suff_len; ++j) {
                    get_word_suflen(gsa[j], nwordi, nsuff_len);
                    if (nsuff_len != suff_len) die("something went wrong!");
                    c = gsa[j]-1 ? dict[gsa[j]-1] : 0;
                    chars.push_back(c);
                    words.push_back(nwordi);
                    same_char = same_char ? (c == pc) : 0;
                    pc = c;
                } // everything seemingly good up till here.
                if ((!any_sa && same_char) || (any_sa && (words.size() == 1)) ) {
                    // print c to bwt after getting all the lengths
                    for (auto word: words)  {
                        for (auto k: get_word_ilist(word)) {
                            if (any_sa) {
                                UPDATE_SA(chars[0], bwsai[k], Difficulty::EASY2);
                            }
                            /*
                            if (update_da) {
                                UPDATE_DA(sa);
                            }
                            */
                            UPDATE_BWT(chars[0]);
                            ++pos;
                            ++easy_cases;
                        }
                    }
                } else {
                    // TODO: maybe a heap will be better? Like in the original
                    // or it's probably faster to sort ahead of time? IDK, must test
                    std::vector<SuffixT> suffs;
                    for (size_t idx = 0; idx < words.size(); ++idx) {
                        // get ilist of each of these words, make a heap
                        auto word = words[idx];
                        for (auto k: get_word_ilist(word)) {
                            suffs.push_back(SuffixT(chars[idx], k, word));
                        }
                    }
                    std::sort(suffs.begin(), suffs.end());
                    for (auto s: suffs) {
                        if (any_sa) {
                            UPDATE_SA(s.bwtc, bwsai[s.bwtp], Difficulty::HARD);
                        }
                        /*
                        if (update_da) {
                            UPDATE_DA(sa);
                        }
                        */
                        UPDATE_BWT(s.bwtc);
                        ++pos;
                        ++hard_cases;
                    }
                }
                chars.clear();
                words.clear();
                next = j;
            }
        }
        if (build_rssa)
            sa_fn(sa_fn_arg(pos-1, psa, RunType::END, Difficulty::EASY1));
        return;
    }

    // void generate_bwt_fm() {
    //     return;
    // }

    private:

    /* run gSACAK on d
     * populates sa, lcp, and dict_idx;
     * this is where the bulk of the algorithm takes its time
     */
    void sort_dict_suffixes(bool build_lcp = true) {
        if (dsize < 1) die("error: dictionary not loaded\n");
        gsa.init_file(fname + "." + EXTGSA, dsize);
        glcp.init_file(fname + "." + EXTGLCP, dsize);
        if (build_lcp)
            gsacak(&dict[0], &gsa[0], &glcp[0], NULL, dsize);
        else { // for when memory is low
            gsacak(&dict[0], &gsa[0], NULL, NULL, dsize);
            // TODO: build FM index over dict
        }
        // make index of dict end positions
        dict_idx = sdsl::bit_vector(dsize, 0);
        for (size_t i = 1; i < dwords + 1; ++i) {
            dict_idx[gsa[i]] = 1;
        }
        dict_idx.init_rs();
    }


    void load_ilist_idx(std::string fname) {
        ReadConType<UIntType> occs(fname + "." + EXTOCC);
        dwords = occs.size();
        uint64_t total_occs = 0;
        for (size_t i = 0; i < occs.size(); ++i)
            total_occs += occs[i];
        ilist_idx = sdsl::bit_vector(total_occs+occs[dwords-1], 0); // TODO: do an assert here
        size_t o = 0;
        for (size_t i = 0; i < occs.size(); ++i) {
            o += occs[i];
            ilist_idx[o-1] = 1;
        }
        ilist_idx.init_rs();
    }

    size_t get_ilist_size(size_t wordi) const {
        auto startpos = wordi ? ilist_idx.select(wordi) + 1 : 0;
        auto endpos = wordi >= dwords ? ilist.size()-1 : ilist_idx.select(wordi+1);
        return endpos - startpos + 1;
    }

    std::vector<size_t> get_word_ilist(size_t wordi) const {
        // get to the end of the previous word's list, then add one to get
        // to the start of the current word
        auto startpos = wordi ? ilist_idx.select(wordi) + 1 : 0;
        auto endpos = wordi >= dwords ? ilist.size()-1 : ilist_idx.select(wordi+1);
        std::vector<size_t> v;
        v.reserve(endpos - startpos + 1);
        for (size_t j = startpos+1; j < endpos+2; ++j)
            v.push_back(ilist[j]);
        return v;
    }

    std::string fname; // prefix fname for storing and loading relevant files
    size_t w=10; // word size of parser
    bool mmapped = false;
    uint64_t dsize; // number of characters in dict
    uint64_t dwords; // number of words in dict
    WriteConType<uint8_t> dict; // dict word array (word ends represented by EndOfWord)
    ReadConType<uint8_t> bwlast; // parse-bwt char associated w/ ilist
    ReadConType<UIntType> ilist; // bwlast positions of dict words
    ReadConType<UIntType> bwsai; // TODO: this might need a separate UIntType
    WriteConType<UIntType> gsa; // gSA of dict words
    WriteConType<IntType> glcp; // gLCP of dict words
    bv_rs<> ilist_idx; // bitvec w/ 1 on ends of dict word occs in ilist
    bv_rs<> dict_idx; // bitvec w/ 1 on word end positions in dict
    bool build_sa = false;
    bool build_rssa = false;
    bool any_sa = false;
    bool verbose = false;
};
}; // namespace end
#endif
