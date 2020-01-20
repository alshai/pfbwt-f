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
    SuffixT(uint8_t c, uint32_t i) : bwtc(c), bwtp(i) {}
    uint8_t bwtc;
    size_t bwtp;
    bool operator<(SuffixT& r);
};

bool SuffixT::operator<(SuffixT& r) {
    return this->bwtp < r.bwtp;
}

struct BwtT {
    BwtT(char ch, uint_t sa, bool b) : c(ch), s(sa), sa(b) {}
    char c; // BWT[i]
    uint_t s; // SA[i]
    bool sa; // whether s should be read or not
};

template<typename IntType, 
         template <typename, typename...> typename ReadConType,
         template <typename, typename...> typename WriteConType
         >
class PrefixFreeBWT {

    public:

    PrefixFreeBWT(std::string prefix, size_t win_size, bool sa = false) :
        fname(prefix),
        w ( win_size),
        dict ( WriteConType<uint8_t>(prefix + "." + EXTDICT)),
        bwlast ( ReadConType<uint8_t>(prefix + "." + EXTBWLST)),
        ilist ( ReadConType<IntType>(prefix + "." + EXTILIST))
    {
        dsize = dict.size();
        load_ilist_idx(prefix);
        if (sa) bwsai = ReadConType<IntType>(prefix + "." + EXTBWSAI);
    }

#define get_word_suflen(i, d, s) \
    d = dict_idx.rank(i); \
    s = d>=dwords ? dsize-i : dict_idx.select(d+1) - i;

     /* uses LCP of dict to build BWT (less memory, more time)
     */
    template<typename F>
    void generate_bwt_lcp(F bwt_processor, bool build_sa=false) {
        sort_dict_suffixes(true); // build gSA and gLCP of dict
        // start from SA item that's not EndOfWord or EndOfDict
        size_t next, suff_len, wordi, ilist_pos;
        uint8_t bwtc, pbwtc = 0;
        size_t easy_cases = 0, hard_cases = 0;
        for (size_t i = dwords+w+1; i<dsize; i=next) {
            next = i+1;
            // fprintf(stdout, "%lu\n", sa[i]);
            get_word_suflen(gsa[i], wordi, suff_len);
            if (suff_len <= w) continue; // ignore small suffixes
            // full word case
            if (gsa[i] == 0 || dict_idx[gsa[i]-1] == 1) {
                auto word_ilist = get_word_ilist(wordi);
                for (auto j: word_ilist) {
                    bwtc = bwlast[j];
                    // TODO: do SA/DA/etc stuff here,
                    if (build_sa) {
                        // get SA
                        bwt_processor(BwtT(bwtc, 0, false));
                    }  else {
                        bwt_processor(BwtT(bwtc, 0, false));
                    }
                    pbwtc = bwtc; // for SA
                }
                ++easy_cases;
            } else { // hard case!
                // look at all the sufs that share LCP[suf]==this_suffixlen
                size_t nwordi, nsuff_len;
                std::vector<uint8_t> chars;
                std::vector<uint64_t> words;
                uint8_t c, pc = dict[gsa[i]-1];
                chars.push_back(pc);
                words.push_back(wordi);
                bool same_char = true;
                size_t j;
                for (j = i + 1; j < dsize && glcp[j] >= (int_t) suff_len; ++j) {
                    get_word_suflen(gsa[j], nwordi, nsuff_len);
                    if (nsuff_len != suff_len) die("something went wrong!");
                    c = dict[gsa[j]-1];
                    chars.push_back(c);
                    words.push_back(nwordi);
                    same_char = same_char ? (c == pc) : 0;
                    pc = c;
                } // everything seemingly good up till here.
                next = j;
                if (same_char || (build_sa && (words.size() == 1)) ) {
                    // print c to bwt after getting all the lengths
                    for (auto word: words)  {
                        size_t nilist = get_ilist_size(word);
                        for (size_t k = 0; k < nilist; ++k) {
                            if (build_sa) {
                                // get SA
                                bwt_processor(BwtT(chars[0], 0, false));
                            } else {
                                bwt_processor(BwtT(chars[0], 0, false));
                            }
                        }
                    }
                    ++easy_cases;
                } else {
                    // TODO: maybe a heap will be better? Like in the original
                    // or it's probably faster to sort ahead of time? IDK, must test
                    std::vector<SuffixT> suffs;
                    for (size_t idx = 0; idx < words.size(); ++idx) {
                        // get ilist of each of these words, make a heap
                        for (auto word: get_word_ilist(words[idx])) {
                            suffs.push_back(SuffixT(chars[idx], word));
                        }
                    }
                    std::sort(suffs.begin(), suffs.end());
                    for (auto s: suffs) {
                        if (build_sa) {
                            // get SA
                            bwt_processor(BwtT(s.bwtc, 0, false));
                        } else {
                            bwt_processor(BwtT(s.bwtc, 0, false));
                        }
                    }
                    ++hard_cases;
                }
                chars.clear();
                words.clear();
            }
        }
        return;
    }

    void generate_bwt_fm() {
        return;
    }

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
        dict[0] = 0; // TODO: I would rather not write *anything* to dict
    }


    void load_ilist_idx(std::string fname) {
        ReadConType<IntType> occs(fname + "." + EXTOCC);
        dwords = occs.size();
        int total_occs = 0;
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

    size_t get_ilist_size(size_t wordi) {
        auto startpos = wordi ? ilist_idx.select(wordi) + 1 : 0;
        auto endpos = wordi >= dwords ? ilist.size()-1 : ilist_idx.select(wordi+1);
        return endpos - startpos + 1;
    }

    std::vector<size_t> get_word_ilist(size_t wordi) {
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
    ReadConType<IntType> ilist; // bwlast positions of dict words
    ReadConType<IntType> bwsai; // TODO: this might need a separate IntType
    WriteConType<uint_t> gsa; // gSA of dict words
    WriteConType<int_t> glcp; // gLCP of dict words
    bv_rs<> ilist_idx; // 1 on ends of dict word occs in ilist
    bv_rs<> dict_idx; // 1 on word end positions in dict
};
}; // namespace end
#endif
