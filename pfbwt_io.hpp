#include <cstdio>
#include <cinttypes>
#include <vector>
#include <string>
#include <tuple>
#include "file_wrappers.hpp"
#include "parse-f.hpp"
extern "C" {
#ifndef AC_KSEQ_H
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
}

size_t get_fasta_length(std::string fname) {
    FILE* fp = fopen(fname.data(), "r");
    char* line = NULL;
    size_t length = 0, n;
    int k;
    while ((k = getline(&line, &n, fp)) >= 0) {
        if (line[0] != '>') {
            length += line[k-1]=='\n' ? k - 1 : k;
        }
    }
    fclose(fp);
    return length;
}

std::vector<std::pair<std::string, size_t>> get_fasta_lengths(std::string fname) {
    std::vector<std::pair<std::string, size_t>> v;
    gzFile fp = gzopen(fname.data(), "r");
    if (fp == NULL) die("failed to open file!\n");
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        v.push_back({seq->name.s, seq->seq.l});
    }
    return v;
}

template<typename T>
std::vector<T> vec_from_file(std::string path) {
    std::vector<T> vec;
    size_t size = get_file_size(path.data());
    size_t nelems = size / sizeof(T);
    FILE* fp = fopen(path.data(), "rb");
    vec.resize(nelems);
    if (fread(&(vec.data())[0], sizeof(T), vec.size(), fp) != nelems) {
        exit(1);
    }
    fclose(fp);
    return vec;
}

// assuming dict is appended by EndOfDict, delimited by EndOfWord
std::vector<std::string> load_dict(std::string dict_fname) {
    FILE* dfp = std::fopen(dict_fname.data(), "rb");
    if (dfp == NULL) die("bad dict file\n");
    int c;
    std::vector<std::string> dvec;
    std::string word;
    while((c = fgetc(dfp)) != EOF) {
        if (c == EndOfDict) {
            break;
        } else if (c == EndOfWord) {
            dvec.push_back(word);
            word.clear();
        } else {
            word.append(1, c);
        }
    }
    fclose(dfp);
    return dvec;
}

template<typename UIntType>
std::vector<UIntType> load_occs(std::string occ_fname) {
    FILE* ofp = fopen(occ_fname.data(), "rb");
    if (ofp == NULL) die("bad occs file\n");
    UIntType occ;
    std::vector<UIntType> occs;
    while (fread(&occ, sizeof(UIntType), 1, ofp) == 1) {
        occs.push_back(occ);
    }
    fclose(ofp);
    return occs;
}

template<typename IntType>
std::vector<IntType> load_parse_ranks(std::string parse_fname) {
    FILE* pfp = fopen(parse_fname.data(), "rb");
    if (pfp == NULL) die("bad parse file\n");
    IntType prank;
    std::vector<IntType> parse_ranks;
    while (fread(&prank, sizeof(IntType), 1, pfp) == 1) {
        parse_ranks.push_back(prank);
    }
    fclose(pfp);
    return parse_ranks;
}

template<typename UIntType>
std::vector<UIntType> load_sai(std::string sai_fname) {
    return vec_from_file<UIntType>(sai_fname);
}

std::vector<char> load_last(std::string last_fname) {
    return vec_from_file<char>(last_fname);
}

std::vector<const char*> load_parse(std::string parse_fname) {
    std::vector<const char*> parse;
    FILE* fp = fopen(parse_fname.data(), "r");
    char* line = NULL;
    size_t n;
    int k;
    while ((k = getline(&line, &n, fp)) >= 0) {
        if (line[k-1] == '\n') line[k-1] = '\0';
        const char* l = strdup(line);
        if (l == NULL) exit(1);
        parse.push_back(l);
    }
    return parse;
}

pfbwtf::Parser<> load_parser(std::string prefix, pfbwtf::ParserParams p) {
    // using UIntType = pfbwtf::Parser<>::UIntType;
    using IntType = pfbwtf::Parser<>::IntType;
    auto dict = load_dict(prefix + ".dict");
    auto parse_ranks = load_parse_ranks<IntType>(prefix + ".parse");
    return pfbwtf::Parser<>(p, dict, std::move(parse_ranks));
}

pfbwtf::Parser<> parse_from_fasta(std::string fasta_fname, pfbwtf::ParserParams p) {
    pfbwtf::Parser<> parser(p);
    parser.add_fasta(fasta_fname);
    parser.finalize();
    return parser;
}
