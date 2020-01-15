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
#include "kseq.h"
extern "C" {
KSEQ_INIT(gzFile, gzread);
#include "utils.h"
#include "gsa/gsacak.h"
}


namespace pfbwtf {

struct ParseArgs {
    std::string in_fname;
    size_t w = 10;
    size_t p = 100;
    bool sai = false;
};

template<typename T>
struct Freq {
    Freq(T x) : n(x) {}
    T n = 0;
    T r = 0;
};


template <typename IntType, typename Hasher>
struct Parser {

    using FreqMap = std::map<std::string, Freq<IntType>>;

    public:

    Parser(size_t wsize, size_t pmod) {
        set_w(wsize);
        set_p(pmod);
    }

    void parse_fasta(const char* fname, const bool sai = false) {
        FILE* sa_fp = NULL;
        if (sai) sa_fp = open_aux_file(fname, EXTSAI, "wb");
        // FILE* old_parse_fp = open_aux_file(fname, "parse2", "w");
        gzFile fp = gzopen(fname, "r");
        kseq_t* seq = kseq_init(fp);
        int l;
        size_t nseqs(0);
        IntType pos = 0;
        std::string phrase;
        phrase.append(1, Dollar);
        Hasher hf(w);
        while (( l = kseq_read(seq) ) >= 0) {
            for (size_t i = 0; i < seq->seq.l; ++i) {
                phrase.append(1, seq->seq.s[i]);
                hf.update(seq->seq.s[i]);
                if (hf.hashvalue() % p == 0) {
                    // fprintf(old_parse_fp, "%s\n", phrase.data());
                    pos = pos ? pos+phrase.size()-w : phrase.size()-1;
                    process_phrase(phrase, sa_fp, &pos);
                    phrase.erase(0, phrase.size()-w);
                }
            }
            ++nseqs;
        }
        phrase.append(w, Dollar);
        // fprintf(old_parse_fp, "%s\n", phrase.data());
        pos = pos ? pos+phrase.size()-w : phrase.size()-1;
        process_phrase(phrase, sa_fp, &pos);

        kseq_destroy(seq);
        gzclose(fp);
        // if (fclose(old_parse_fp)) die("error closing PARSE2 file\n");
        // else fprintf(stderr, "PARSE2 file written to %s.parse2\n", fname);
        if (sai) {
            if (fclose(sa_fp)) die("error closing SAI file\n");
            else (fprintf(stderr, "SAI file written to %s.%s\n", fname, EXTSAI));
        }
        return;
    }

    // assigns lexicographic rankings to items in dictionary
    // and creates occ array
    void update_dict(const char* fname = NULL) {
        std::vector<const char*> dict_phrases;
        dict_phrases.reserve(freqs.size());
        for (auto it = freqs.begin(); it != freqs.end(); ++it) {
            dict_phrases.push_back(it->first.data());
        }
        fprintf(stderr, "sorting dict\n");
        std::sort(dict_phrases.begin(), dict_phrases.end(),
                [](const char* l, const char* r) { return strcmp(l, r) <= 0; });
        fprintf(stderr, "writing dict and occ file...\n");
        FILE* dict_fp = NULL;
        FILE* occ_fp = NULL;
        if (fname != NULL) {
            dict_fp = open_aux_file(fname, EXTDICT, "wb");
            occ_fp  = open_aux_file(fname, EXTOCC, "wb");
        }
        occs.clear();
        occs.reserve(dict_phrases.size());
        size_t rank = 1;
        for (auto x: dict_phrases) {
            // access dictionary and write occ to occ file
            auto& wf = freqs.at(x);
            wf.r = rank++;
            occs.push_back(wf.n);
            if (fname != NULL) {
                if (fwrite(x, 1, strlen(x), dict_fp) != strlen(x))
                    die("Error writing to DICT file\n");
                if (fputc(EndOfWord, dict_fp) == EOF)
                    die("Error writing EndOfWord to DICT file");
                if (fwrite(&wf.n, sizeof(wf.n), 1, occ_fp) != 1)
                    die("Error writing to OCC file\n");
            }
        }
        if (fname != NULL) {
            if (fputc(EndOfDict, dict_fp) == EOF) die("Error writing EndOfDict to DICT file");
            if (fclose(occ_fp)) die("Error closing OCC file");
            else fprintf(stderr, "OCC written to %s.%s\n", fname, EXTOCC);
            if (fclose(dict_fp)) die("Error closing DICT file");
            else fprintf(stderr, "DICT written to %s.%s\n", fname, EXTDICT);
        }
    }

    void generate_parse_ranks() {
        if (!freqs.size()) {
            fprintf(stderr, "dictionary not created yet! run update_dict\n");
            exit(1);
        }
        clear_parse_ranks();
        parse_ranks.reserve(parse.size());
        occs.reserve(parse.size());
        for (auto phrase: parse) {
            auto wf = freqs.at(std::string(phrase));
            parse_ranks.push_back(wf.r);
        }
    }

    /* dumps lexicographic ranks of phrases in the parse
     * to fname.parse
     */
    void dump_parse(const char* fname) {
        if (!parse_ranks.size()) {
            generate_parse_ranks();
        }
        size_t n = parse_ranks.size();
        FILE* parse_fp = open_aux_file(fname, EXTPARSE, "wb");
        n = parse_ranks[n-1] ? n : n-1;
        if (fwrite(parse_ranks.data(), sizeof(parse_ranks[0]), n, parse_fp) != n)
            die("error writing to PARSE");
        if (fclose(parse_fp)) die("Error closing PARSE file");
        else fprintf(stderr, "PARSE written to %s.%s\n", fname, EXTPARSE);
    }

    void dump_last(const char* fname) {
        FILE* last_fp = open_aux_file(fname, EXTLST, "wb");
        if (fwrite(last.data(), sizeof(last[0]), last.size(), last_fp) != last.size())
            die("error writing to LAST");
        if (fclose(last_fp)) die("error closing LAST file\n");
        else fprintf(stderr, "LAST file written to %s.%s\n", fname, EXTLST);
    }

    void dump_occs(const char* fname) {
        fprintf(stderr, "writing occ file...\n");
        FILE* occ_fp = open_aux_file(fname, EXTOCC, "wb");
        if (fwrite(occs.data(), sizeof(occs[0]), occs.size(), occ_fp) != occs.size())
            die("error writing to OCC");
        if (fclose(occ_fp)) die("Error closing OCC file");
        else fprintf(stderr, "OCC written to %s.%s\n", fname, EXTOCC);
    }

    void clear_dict() {
        freqs.clear();
    }

    void clear_parse() {
        parse.clear();
    }

    void clear_parse_ranks() {
        parse_ranks.clear();
    }

    void clear_occ() {
        occs.clear();
    }

    void clear_last() {
        last.clear();
    }

    void clear() {
        clear_dict();
        clear_parse();
        clear_parse_ranks();
        clear_last();
    }

    size_t set_w(size_t x) {
        if (x < 32) w = x;
        else {
            fprintf(stderr, "window size w must be < 32!\n");
            exit(1);
        }
        fprintf(stderr, "w reset. clearing Parser\n");
        clear();
        return w;
    }

    size_t set_p(size_t x) {
        p = x;
        fprintf(stderr, "p reset. clearing Parser\n");
        clear();
        return p;
    }

    // generates bwlast and ilist
    void bwt_of_parse(const char* fname=NULL, bool get_sai=false) {
        size_t n; // size of parse_ranks, minus the last EOS character
        // TODO: add option for using mmap, or turn into sdsl::int_vector
        std::vector<IntType> sai; 
        FILE* sa_fp;
        bwsai.clear();
        // TODO: support large parse sizes 
        if (parse_ranks.size() > 0x7FFFFFFE) {
            fprintf(stderr, "currently no support for texts w/ > 2^31-2 phrases\n");
            exit(1);
        }
        // add EOS if it's not already there
        if (parse_ranks[parse_ranks.size()-1]) {
            n = parse_ranks.size();
            parse_ranks.push_back(0);
        } else n = parse_ranks.size() - 1;
        // TODO: calculate k
        size_t k = 0;
        for (size_t i = 0; i < n; ++i) {
            k = parse_ranks[i] > k ? parse_ranks[i] : k;
        }
        if (get_sai) {
            // TODO: use sdsl int-vector instead
            if (!fname) die("no filename specified\n");
            sa_fp = open_aux_file(fname, EXTSAI, "rb");
            sai.reserve(n);
            fread(&sai[0], sizeof(IntType), n, sa_fp);
            bwsai.reserve(n);
            fclose(sa_fp);
        }
        // compute S.A.
        // we assign instead of reserve in order to be able to use .size()
        assert(sizeof(IntType) >= sizeof(uint_t));
        std::vector<IntType> SA(n+1, 0);
        fprintf(stderr, "Computing S.A. of size %ld over an alphabet of size %ld\n",n+1,k+1);
        int depth = sacak_int(parse_ranks.data(), SA.data(), n+1, k+1);
        if(depth>=0) fprintf(stderr, "S.A. computed with depth: %d\n", depth);
        else die("Error computing the S.A.");
        // transform S.A. to BWT in place
        assert(SA[0] == n);
        SA[0] = parse_ranks[n-1];
        bwlast.clear();
        bwlast.reserve(n+1);
        bwlast.push_back(last[n-2]);
        // TODO: fill bwtsai here
        // if (get_sai) bwsai.push_back(get_myint(&sai[0], n, n-1));
        if (get_sai) bwsai.push_back(sai[n-1]);
        for (size_t i = 1; i < n+1; ++i) {
            if (!SA[i]) {
                SA[i] = 0;
                bwlast.push_back(0);
                // TODO: write to bwtsai
                if (get_sai) bwsai.push_back(0);
            } else {
                if (SA[i] == 1) {
                    bwlast.push_back(last[n-1]);
                } else {
                    bwlast.push_back(last[SA[i]-2]);
                }
                // TODO: write to bwtsai
                // if (get_sai) bwsai.push_back(get_myint(&sai[0], n, SA[i]-1));
                if (get_sai) bwsai.push_back(sai[SA[i]-1]);
                SA[i] = parse_ranks[SA[i] - 1];
            }
        }
        // in bwtparse, S.A. is copied to BWT (ie. text);
        // here, we keep on using S.A., since nothing is begin modified 
        // BWT = BWTSA = S.A.
        std::vector<IntType> F(occs.size()+1, 0);
        F[1] = 1;
        for (size_t i = 2; i < occs.size() + 1; ++i) {
            F[i] = F[i-1] + occs[i-2];
        }
        assert(F[occs.size()] + occs[occs.size()-1] == n+1);
        ilist.resize(n+1, 0);
        for (size_t i = 0; i < n + 1; ++i) {
            ilist[F[SA[i]]++] = i;
        }
        assert(ilist[0]==1);
        assert(SA[ilist[0]] == 0);
    }

    void dump_bwlast(const char* fname) {
        FILE* bwlast_fp = open_aux_file(fname, EXTBWLST, "wb");
        if (fwrite(bwlast.data(), sizeof(char), bwlast.size(), bwlast_fp) != bwlast.size() ) {
            die("could not write BWTLST");
        }
        if (fclose(bwlast_fp)) die("failed writing to BWTLST");
        else fprintf(stderr, "bwlast written to %s.%s\n", fname, EXTBWLST);

    }

    void dump_ilist(const char* fname) {
        FILE* ilist_fp = open_aux_file(fname, EXTILIST, "wb");
        if (fwrite(ilist.data(), sizeof(ilist[0]), ilist.size(), ilist_fp) != ilist.size() ) {
            die("could not write ILIST");
        }
        if (fclose(ilist_fp)) die("failed writing to EXTILIST");
        else fprintf(stderr, "ilist written to %s.%s\n", fname, EXTILIST);
    }

    void dump_bwsai(const char* fname) {
        FILE* sa_fp = open_aux_file(fname, EXTBWSAI, "wb");
        if (fwrite(bwsai.data(), 1, bwsai.size(), sa_fp) != bwsai.size()) { 
            die("could not write BWSAI");
        }
        if (fclose(sa_fp)) die("failed writing to BWSAI");
        else fprintf(stderr, "bwsai written to %s.%s.", fname, EXTBWSAI);
    }

    private:

    void inline process_phrase(const std::string& phrase, FILE* sa_fp=NULL, IntType* pos=NULL) {
        auto ret = freqs.insert({phrase, Freq<IntType>(1)});
        if (!ret.second) ret.first->second.n += 1;
        parse.push_back(ret.first->first.data());
        last.push_back(phrase[phrase.size()-w-1]);
        if (sa_fp != NULL) {
            if(fwrite(pos, sizeof(IntType), 1, sa_fp) != 1)
                die("error writing to .sai file\n");
        }
    }

    FreqMap freqs;
    std::vector<const char*> parse;
    std::vector<IntType> occs;
    std::vector<IntType> parse_ranks;
    std::vector<char> last;
    std::vector<char> bwlast;
    // sdsl::int_vector< bwlast;
    std::vector<IntType> ilist;
    std::vector<IntType> bwsai;
    size_t w, p;
};
}; // namespace end

#endif
