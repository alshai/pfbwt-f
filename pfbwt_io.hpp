#include <cstdio>
#include <cinttypes>
#include <vector>
#include <string>
#include "file_wrappers.hpp"
#include "parse-f.hpp"

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

pfbwtf::Parser<> load_parser(std::string prefix) {
    using UIntType = pfbwtf::Parser<>::UIntType;
    using IntType = pfbwtf::Parser<>::IntType;
    pfbwtf::ParserParams p;
    std::vector<std::string> dict = load_dict(prefix + ".dict");
    auto parse_ranks = load_parse_ranks<IntType>(prefix + ".parse");
    auto occs = load_occs<UIntType>(prefix + ".occ");
    return pfbwtf::Parser<>(p, dict, occs, parse_ranks);
}

pfbwtf::Parser<> parse_from_fasta(std::string fasta_fname) {
    pfbwtf::ParserParams p;
    pfbwtf::Parser<> parser(p);
    parser.add_fasta(fasta_fname);
    parser.finalize();
    return parser;
}
