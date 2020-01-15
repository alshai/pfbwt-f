#ifndef PFBWT_HPP
#define PFBWT_HPP

#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>
#include <fcntl.h>
#include "sdsl/bit_vectors.hpp"
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

struct BwtT {
    BwtT(char ch, uint_t sa, bool b) : c(ch), s(sa), sa(b) {}
    char c; // BWT[i]
    uint_t s; // SA[i]
    bool sa; // whether s should be read or not
};

bool SuffixT::operator<(SuffixT& r) {
    return this->bwtp < r.bwtp;
}

template<typename IntType, 
         template <typename, typename...> typename ReadConType,
         template <typename, typename...> typename WriteConType
         >
class PrefixFreeBWT {

    public:

    PrefixFreeBWT(std::string prefix, size_t win_size, bool sa = false) :
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
    d = dict_idx_rank(i); \
    s = d>=dwords ? dsize-i : dict_idx_select(d+1) - i;

     /* uses LCP of dict to build BWT (less memory, more time)
     */
    template<typename F>
    void generate_bwt_lcp(F bwt_processor, bool build_sa=false) {
        sort_dict_suffixes(true); // build SA and lcp
        // start from SA item that's not EndOfWord or EndOfDict
        size_t next, suff_len, wordi, ilist_pos;
        uint8_t bwtc, pbwtc = 0;
        size_t easy_cases = 0, hard_cases = 0;
        for (size_t i = dwords+w+1; i<dsize; i=next) {
            next = i+1;
            // fprintf(stdout, "%lu\n", sa[i]);
            get_word_suflen(sa[i], wordi, suff_len);
            if (suff_len <= w) continue; // ignore small suffixes
            // full word case
            if (sa[i] == 0 || dict_idx[sa[i]-1] == 1) {
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
                uint8_t c, pc = dict[sa[i]-1];
                chars.push_back(pc);
                words.push_back(wordi);
                bool same_char = true;
                size_t j;
                for (j = i + 1; j < dsize && lcp[j] >= (int_t) suff_len; ++j) {
                    get_word_suflen(sa[j], nwordi, nsuff_len);
                    if (nsuff_len != suff_len) die("something went wrong!");
                    c = dict[sa[j]-1];
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
        sa.resize(dsize);
        lcp.resize(dsize);
        if (build_lcp)
            gsacak(&dict[0], &sa[0], &lcp[0], NULL, dsize);
        else { // for when memory is low
            gsacak(&dict[0], &sa[0], NULL, NULL, dsize);
            // TODO: build FM index over dict
        }
        // make index of dict end positions
        dict_idx = sdsl::bit_vector(dsize, 0);
        for (size_t i = 1; i < dwords + 1; ++i) {
            dict_idx[sa[i]] = 1;
        }
        sdsl::util::init_support(dict_idx_rank, &dict_idx);
        sdsl::util::init_support(dict_idx_select, &dict_idx);
        dict[0] = 0;
    }

    /*
    void load_files(std::string prefix, bool sa=false) {
        dict = ConType<uint8_t>(prefix + "." + EXTDICT);
        bwast = ConType<uint8_t>(prefix + "." + EXTBWLST);
        ilist = ConType<IntType>(prefix + "." + EXTILIST);
        load_ilist_idx(prefix);
        if (sa) bwsai = ConType<IntType>(prefix + "." + EXTBWSAI)
    }
    */

    void load_ilist_idx(std::string fname) {
        ReadConType<IntType> occs(fname + "." + EXTOCC);
        /*
        FILE* occ_fp = open_aux_file(fname.data(), EXTOCC, "rb");
        // get occs size
        fseek(occ_fp, 0, SEEK_END);
        size_t osize = ftell(occ_fp);
        if (osize % 4 != 0) die("invalid occ file");
        dwords = osize / sizeof(IntType);
        rewind(occ_fp);
        std::vector<uint32_t> occs(dwords);
        if (fread(&occs[0], sizeof(IntType), occs.size(), occ_fp) != dwords) {
            die("error reading occs");
        }
        int total_occs = 0;
        for (auto o: occs) total_occs += o;
        */
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
        sdsl::util::init_support(ilist_idx_rank,   &ilist_idx);
        sdsl::util::init_support(ilist_idx_select, &ilist_idx);
        // fclose(occ_fp);
    }

    size_t get_ilist_size(size_t wordi) {
        auto startpos = wordi ? ilist_idx_select(wordi) + 1 : 0;
        auto endpos = wordi >= dwords ? ilist.size()-1 : ilist_idx_select(wordi+1);
        return endpos - startpos + 1;
    }

    std::vector<size_t> get_word_ilist(size_t wordi) {
        // get to the end of the previous word's list, then add one to get
        // to the start of the current word
        auto startpos = wordi ? ilist_idx_select(wordi) + 1 : 0;
        auto endpos = wordi >= dwords ? ilist.size()-1 : ilist_idx_select(wordi+1);
        std::vector<size_t> v;
        v.reserve(endpos - startpos + 1);
        for (size_t j = startpos+1; j < endpos+2; ++j)
            v.push_back(ilist[j]);
        // fprintf(stderr, "%lu %lu\n", startpos+1, endpos+1);
        return v;
    }

    size_t w=10; // word size of parser
    bool mmapped = false;
    uint64_t dsize; // number of characters in dict
    uint64_t dwords; // number of words in dict
    WriteConType<uint8_t> dict; // dict word array (word ends represented by EndOfWord)
    ReadConType<uint8_t> bwlast; // parse-bwt char associated w/ ilist
    ReadConType<IntType> ilist; // bwlast positions of dict words
    ReadConType<IntType> bwsai; // TODO: this might need a separate IntType
    sdsl::bit_vector ilist_idx; // 1 on ends of dict word occs in ilist
    sdsl::bit_vector dict_idx; // 1 on word end positions in dict
    sdsl::bit_vector::rank_1_type   ilist_idx_rank;
    sdsl::bit_vector::rank_1_type   dict_idx_rank;
    sdsl::bit_vector::select_1_type ilist_idx_select;
    sdsl::bit_vector::select_1_type dict_idx_select;
    // TODO: implement read-writable ConType so we can handle these
    //       this is important bc memory bottleneck is here
    std::vector<uint_t> sa; // gSA of dict words
    std::vector<int_t> lcp; // gLCP of dict words
};
}; // namespace end
#endif
