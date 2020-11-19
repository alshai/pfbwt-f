#include <cstdio>
#include <cinttypes>
#include <vector>
#include <string>
#include <tuple>
#include <sys/stat.h>
#include "file_wrappers.hpp"
#include "pfparser.hpp"
extern "C" {
#include "utils.h"
#ifndef AC_KSEQ_H
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
}

namespace pfbwtf {

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
void vec_to_file(const std::vector<T>& vec, std::string fname) {
    FILE* fp = fopen(fname.data(), "wb");
    if (fwrite(vec.data(), sizeof(T), vec.size(), fp) != vec.size() ) {
        die("could not write file");
    }
    fclose(fp);
}

template<>
void vec_to_file<std::string>(const std::vector<std::string>& vec, std::string fname) {
    FILE* fp = fopen(fname.data(), "w");
    for (auto v: vec) {
        fprintf(fp, "%s\n", v.data());
    }
    fclose(fp);
}

template<typename T>
void vec_to_file(const std::vector<T>& vec, size_t nelems, std::string fname) {
    FILE* fp = fopen(fname.data(), "wb");
    if (fwrite(vec.data(), sizeof(T), nelems, fp) != nelems ) {
        die("could not write file");
    }
    fclose(fp);
}

void dict_to_file(const std::vector<const char*>& phrases, std::string fname) {
    FILE* dict_fp = fopen(fname.data(), "wb");
    if (dict_fp == NULL) die("unable to open dict file");
    for (auto phrase: phrases) {
        if (fwrite(phrase, 1, strlen(phrase), dict_fp) != strlen(phrase))
            die("Error writing to DICT file\n");
        if (fputc(EndOfWord, dict_fp) == EOF)
            die("Error writing EndOfWord to DICT file");
    }
    if (fputc(EndOfDict, dict_fp) == EOF) die("Error writing EndOfDict to DICT file");
    if (fclose(dict_fp)) die("Error closing DICT file");
}

template<typename T>
std::vector<T> vec_from_file(std::string path) {
    std::vector<T> vec;
    size_t size = get_file_size(path.data());
    size_t nelems = size / sizeof(T);
    FILE* fp = fopen(path.data(), "rb");
    if (fp == NULL) {
        fprintf(stderr, "vec_from_file: error opening %s\n", path.data());
        exit(1);
    }
    vec.resize(nelems);
    if (fread(&(vec.data())[0], sizeof(T), vec.size(), fp) != nelems) {
        fprintf(stderr, "vec_from_file: error reading from %s\n", path.data());
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
        if (l == NULL) {
            fprintf(stderr, "error duplicating string\n");
            exit(1);
        }
        parse.push_back(l);
    }
    return parse;
}

template<typename U>
std::pair<std::vector<std::string>, std::vector<U>> load_doc_info(std::string fname) {
    std::vector<std::string> names;
    std::vector<U> starts;
    FILE* fp = fopen(fname.data(), "r");
    if (fp == NULL) {
        fprintf(stderr, "error opening doc file %s", fname.data());
        exit(1);
    }
    char* line = NULL;
    char* tok = NULL, *end = NULL;
    size_t n = 0;
    int l;
    const char* delim = " ";
    char* saveptr;
    while ((l = getline(&line, &n, fp) >= 0)) {
        tok = strtok_r(line, delim, &saveptr);
        if (tok == NULL) {
            fprintf(stderr, "error parsing name in doc file line\n");
            exit(1);
        }
        names.push_back(tok);
        tok = strtok_r(NULL, delim, &saveptr);
        if (tok == NULL) {
            fprintf(stderr, "error parsing size in doc file line: %s\n", line);
            exit(1);
        }
        starts.push_back(std::strtoul(tok, &end, 10));
    }
    free(line);
    return std::make_pair(names, starts);
}

/* loads parser from .dict and .parse files */
pfbwtf::PfParser<> load_parser(std::string prefix, pfbwtf::PfParserParams p) {
    using UIntType = pfbwtf::PfParser<>::UIntType;
    using IntType = pfbwtf::PfParser<>::IntType;
    auto dict = load_dict(prefix + ".dict");
    auto parse_ranks = load_parse_ranks<IntType>(prefix + ".parse");
    if (p.store_docs) {
        auto doc_pair = load_doc_info<UIntType>(prefix + ".docs");
        return pfbwtf::PfParser<>(p, dict, std::move(parse_ranks), std::move(doc_pair.second), std::move(doc_pair.first));
    } else {
        return pfbwtf::PfParser<>(p, dict, std::move(parse_ranks));
    }
}

template<typename U>
void docs_to_file(std::string fname, const std::vector<std::string>& doc_names, const std::vector<U>& doc_starts) {
    FILE* doc_fp = fopen(fname.data(), "w");
    for (size_t i = 0; i < doc_starts.size(); ++i) {
        fprintf(doc_fp, "%s %lu\n", doc_names[i].data(), doc_starts[i]);
    }
    fclose(doc_fp);
}

/* saves parser to .dict, .occ, and .parse files (and .docs if applicable)*/
void save_parser(const pfbwtf::PfParser<>& parser, std::string prefix) {
    std::string dict_fname = prefix + ".dict";
    std::string occ_fname = prefix + ".occ";
    std::string n_fname = prefix + ".n";
    std::string parse_ranks_fname = prefix + ".parse";
    dict_to_file(parser.get_sorted_phrases(), dict_fname);
    vec_to_file(parser.get_occs(), occ_fname);
    vec_to_file(parser.get_parse_ranks(), parser.get_parse_size(), parse_ranks_fname);
    if (parser.get_params().store_docs) {
        std::string docs_fname = prefix + ".docs";
        docs_to_file(docs_fname, parser.get_doc_names(), parser.get_doc_starts());
    }
    std::FILE* n_fp = fopen(n_fname.data(), "w");
    fprintf(n_fp, "%lu\n", parser.get_n());
    fclose(n_fp);
}

pfbwtf::PfParser<> parse_from_fasta(std::string fasta_fname, pfbwtf::PfParserParams p) {
    pfbwtf::PfParser<> parser(p);
    parser.add_fasta(fasta_fname);
    parser.finalize();
    return parser;
}

int parse_files_exist(std::string prefix) {
    std::string dict_prefix = prefix + ".dict";
    std::string parse_prefix = prefix + ".parse";
    struct stat buffer;
    return (!stat(dict_prefix.data(), &buffer) &&  !stat(parse_prefix.data(), &buffer));
}

int file_exists(std::string fname) {
    struct stat buffer;
    return (!stat(fname.data(), &buffer));
}

PfParser<> load_or_generate_parser_w_log(std::string prefix, PfParserParams params, FILE* fp = stderr) {
    PfParser<> parser;
    if (parse_files_exist(prefix)) {
        fprintf(fp, "loading %s from file\n", prefix.data());
        parser += load_parser(prefix, params);
    } else {
        fprintf(fp, "generating parse for %s\n", prefix.data());
        if (file_exists(prefix)) {
            parser += parse_from_fasta(prefix, params);
        } else {
            // TODO: figure out how to exit a threaded program here and do cleanup
            fprintf(fp, "ERROR: %s not found, cannot add it to parse!\n", prefix.data());
        }
    }
    return parser;
}

void save_parse_bwt(PfParser<>& parser, std::string output, bool sa = false) {
        parser.bwt_of_parse(
                [&](const std::vector<char>& bwlast,
                    const std::vector<typename PfParser<>::UIntType>& ilist,
                    const std::vector<typename PfParser<>::UIntType>& bwsai)
                {
                    vec_to_file<char>(bwlast, output + ".bwlast");
                    vec_to_file<PfParser<>::UIntType>(ilist, output + ".ilist");
                    if (sa) vec_to_file<typename PfParser<>::UIntType>(bwsai, output + ".bwsai");
                });
}

} // namespace
