#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "file_wrappers.hpp"
#include "pfparser.hpp"
#include "pfbwt_io.hpp"
extern "C" {
#include <utils.h>
}

bool parser_cmp(const pfbwtf::PfParser<>& lhs, const pfbwtf::PfParser<>& rhs, std::string msg, FILE* log) {
    // if (lhs.get_pos() != rhs.get_pos()) {
    //     fprintf(log, "%s: %s: pos mismatch %lu vs %lu\n", msg.data(), __func__, lhs.get_pos(), rhs.get_pos());
    //     return false;
    // }
    if (lhs.get_parse_ranks().size() != rhs.get_parse_ranks().size()) {
        fprintf(log, "%s: %s: parse_ranks_ size mismatch %lu vs %lu\n", msg.data(), __func__, lhs.get_parse_ranks().size(), rhs.get_parse_ranks().size());
        for (auto p: lhs.get_parse_ranks()) { fprintf(log, "%d ", p); } fprintf(log, "\n");
        for (auto p: rhs.get_parse_ranks()) { fprintf(log, "%d ", p); } fprintf(log, "\n");
        for (auto p: lhs.get_parse()) { fprintf(log, "%s ", p); } fprintf(log, "\n");
        for (auto p: rhs.get_parse()) { fprintf(log, "%s ", p); } fprintf(log, "\n");
        return false;
    }
    if (lhs.get_last().size() != rhs.get_last().size())  {
        fprintf(log, "%s: %s: last_ size mismatch %lu vs %lu\n", msg.data(), __func__, lhs.get_last().size(), rhs.get_last().size());
        return false;
    }
    if (lhs.get_freqs().size() != rhs.get_freqs().size()) {
        fprintf(log, "%s: %s: freqs_ size mismatch %lu vs %lu\n", msg.data(), __func__, lhs.get_freqs().size(), rhs.get_freqs().size());
        return false;
    }
    if (lhs.get_params().get_sai) {
        if (lhs.get_sai().size() != rhs.get_sai().size()) {
            fprintf(log, "%s: %s: sai_ size mismatch %lu vs %lu\n", msg.data(), __func__, lhs.get_sai().size(), rhs.get_sai().size());
            return false;
        }
    }
    const auto& lhs_freqs = lhs.get_freqs();
    const auto& rhs_freqs = rhs.get_freqs();
    for (const auto& pair: lhs_freqs) {
        const auto it2 = rhs_freqs.find(pair.first);
        if (it2 == rhs_freqs.end()) {
            fprintf(log, "%s: %s: key mismatch (%s)\n", msg.data(), __func__, pair.first.data());
            fprintf(log, "one: "); for (auto k: lhs_freqs) {fprintf(log, "%s ", k.first.data()); } fprintf(log, "\n");
            fprintf(log, "two: "); for (auto k: rhs_freqs) {fprintf(log, "%s ", k.first.data()); } fprintf(log, "\n");
            return false;
        }
        if (it2->second != pair.second) {
            fprintf(log, "%s: %s: item mismatch (%s) %lu %lu vs %lu %lu\n", msg.data(), __func__, pair.first.data(), pair.second.n, pair.second.r, it2->second.n, it2->second.r);
            return false;
        }
    }
    const auto& lhs_parse_ranks = lhs.get_parse_ranks();
    const auto& rhs_parse_ranks = rhs.get_parse_ranks();
    for (size_t i = 0; i < lhs_parse_ranks.size(); ++i) {
        if (lhs_parse_ranks[i] != rhs_parse_ranks[i]){
            fprintf(log, "%s: %s: parse rank mismatch at %lu: %d vs %d\n", msg.data(), __func__, i, lhs_parse_ranks[i], rhs_parse_ranks[i]);
            return false;
        }
    }
    const auto& lhs_last = lhs.get_last();
    const auto& rhs_last = rhs.get_last();
    for (size_t i = 0; i < lhs_last.size(); ++i) {
        if (lhs_last[i] != rhs_last[i]) {
            fprintf(log, "%s: %s: last mismatch at %lu: %c vs %c\n", msg.data(), __func__, i, lhs_last[i], rhs_last[i]);
            return false;
        }
    }
    if (lhs.get_params().get_sai) {
        const auto& lhs_sai = lhs.get_sai();
        const auto& rhs_sai = rhs.get_sai();
        for (size_t i = 0; i < lhs_sai.size(); ++i) {
            if (lhs_sai[i] != rhs_sai[i]) {
                fprintf(log, "%s: %s: sai mismatch at %lu / %lu: %lu vs %lu\n", msg.data(), __func__, i, lhs_sai.size(), lhs_sai[i], rhs_sai[i]);
                for (auto p: lhs.get_sai()) { fprintf(log, "%lu ", p); } fprintf(log, "\n");
                for (auto p: rhs.get_sai()) { fprintf(log, "%lu ", p); } fprintf(log, "\n");
                return false;
            }
        }
    }
    return true;
}

// makes sure that the loaded file's aux data structures are properly generated
bool parser_test_load(FILE* log) {
    // load parser
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    std::string prefix = "tests/random_examples/random.all.fa";
    pfbwtf::PfParser<> p_loaded(load_parser(prefix, params));
    auto lhs_parse = load_parse(prefix + ".parse_"); // NOTE: need to free elements here
    auto rhs_parse = p_loaded.get_parse();
    if (lhs_parse.size() != rhs_parse.size()) {
        fprintf(log, "%s: parse sizes unequal\n", __func__);
        return false;
    }
    for (size_t i = 0; i < lhs_parse.size(); ++i) {
        if (strcmp(lhs_parse[i], rhs_parse[i])) {
            fprintf(log, "%s: parse unequal at %lu (%s vs %s)\n", __func__, i, lhs_parse[i], rhs_parse[i]);
            return false;
        }
    }
    auto lhs_last = load_last(prefix + ".last");
    auto rhs_last = p_loaded.get_last();
    if (lhs_last.size() != rhs_last.size()) {
        fprintf(log, "%s: last sizes unequal\n", __func__);
        return false;
    }
    for (size_t i = 0; i < lhs_last.size(); ++i)  {
        if (lhs_last[i] != rhs_last[i]) {
            fprintf(log, "%s: last unequal at %lu (%c vs %c)\n", __func__, i, lhs_last[i], rhs_last[i]);
            return false;
        }
    }

    auto lhs_sai = load_sai<uint64_t>(prefix + ".sai");
    auto rhs_sai = p_loaded.get_sai();
    if (lhs_sai.size() != rhs_sai.size()) {
        fprintf(log, "%s: sai sizes unequal\n", __func__);
        return false;
    }
    for (size_t i = 0; i < lhs_sai.size(); ++i)  {
        if (lhs_sai[i] != rhs_sai[i]) {
            fprintf(log, "%s: sai unequal at %lu: (%lu vs %lu)\n", __func__, i, lhs_sai[i], rhs_sai[i]);
            return false;
        }
    }
    // for (auto s: lhs_parse) free(s);
    return true;
}

bool parser_test_eq(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    pfbwtf::PfParser<> truth(load_parser("tests/random_examples/random.1.fa", params));
    pfbwtf::PfParser<> test;
    test = truth;
    return parser_cmp(truth, test, std::string(__func__), log);
}

// mult. seqs in one fasta file
int parser_test_add_fasta1(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    pfbwtf::PfParser<> truth(load_parser("tests/random_examples/random.all.fa", params));
    pfbwtf::PfParser<> test(parse_from_fasta("tests/random_examples/random.all.fa", params));
    return parser_cmp(truth, test, std::string(__func__), log);
}

// checks if loading seqs from sequence of files matches loading all seqs from a single file
bool parser_test_add_fasta2(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    pfbwtf::PfParser<> truth(load_parser("tests/random_examples/random.all.fa", params));
    pfbwtf::PfParser<> test(params);
    for (int i = 1; i < 251; ++i) {
        std::string prefix = "tests/random_examples/random." + std::to_string(i) + ".fa";
        test.add_fasta(prefix);
    }
    test.finalize();
    return parser_cmp(truth, test, std::string(__func__), log);
}

bool parser_test_pluseq(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    pfbwtf::PfParser<> truth(load_parser("tests/random_examples/random.all.fa", params));
    pfbwtf::PfParser<> test(params);
    for (int i = 1; i < 251; ++i) {
        std::string prefix = "tests/random_examples/random." + std::to_string(i) + ".fa";
        pfbwtf::PfParser<> p(parse_from_fasta(prefix, params));
        test += p;
    }
    test.finalize();
    return parser_cmp(truth, test, std::string(__func__), log);
}

bool parser_test_get_n(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    std::string fname = "tests/random_examples/random.all.fa";
    pfbwtf::PfParser<> p(load_parser(fname, params));
    auto fasta_info = get_fasta_lengths(fname);
    size_t true_n = 0;
    for (auto x: fasta_info) { true_n += x.second + params.w; } // add w bc we implicitly add As
    fprintf(log, "%s: true n == %lu, n == %lu\n", __func__, true_n, p.get_n());
    return true_n == p.get_n();
}

constexpr int NMERGE = 10;
bool parser_test_merge(FILE* log) {
    pfbwtf::PfParserParams params;
    params.get_sai = true;
    pfbwtf::PfParser<> to_merge[250 / NMERGE];
    std::string prefix = "tests/random_examples/random.";
    for (int i = 0; i < 250; i += NMERGE) {
        for (int j = i; j < i + NMERGE; ++j) {
            std::string fname = prefix + std::to_string(j+1) + ".fa";
            to_merge[i/NMERGE] += load_parser(fname, params);
        }
        to_merge[i/NMERGE].finalize(); // TODO: I don't understand why, but this finalize step is very necessary
    }
    pfbwtf::PfParser<> merged;
    for (int i = 0; i < 250 / NMERGE; ++i) {
        merged += to_merge[i];
    }
    merged.finalize(); // generate ranks etc
    pfbwtf::PfParser<> truth(load_parser(prefix + "all.fa", params));
    return parser_cmp(truth, merged, std::string(__func__), log);
}

bool print_test(std::string msg, bool b) {
    fprintf(stdout, "%s: %s\n", msg.data(), b ? "success" : "fail");
    return b;
}

int main() {
    FILE* log = fopen("test_parser.log", "w");
    print_test("load_from_disk", parser_test_load(log));
    print_test("=", parser_test_eq(log));
    print_test("add_fasta1_(one fasta)", parser_test_add_fasta1(log));
    print_test("add_fasta2_(mult. fasta)", parser_test_add_fasta2(log));
    print_test("+=", parser_test_pluseq(log));
    print_test("n", parser_test_get_n(log));
    print_test("merge", parser_test_merge(log));

    fclose(log);
    return 0;
}
