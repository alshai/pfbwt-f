#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <getopt.h>
#include <chrono>
#include "pfbwt-f.hpp"
#include "parse-f.hpp"
#include "hash.hpp"
#include "file_wrappers.hpp"
extern "C" {
#include "utils.h"
}


struct Args {
    std::string in_fname;
    size_t w = 10;
    size_t p = 100;
    bool sa = false;
    bool rssa = false;
    bool mmap = false;
    int parse_only = 0;
};

struct Timer {
    using clock = std::chrono::system_clock;
    using sec = std::chrono::duration<double>;

    Timer(std::string m) : msg(m) {
        start = clock::now();
    }

    ~Timer() {
        sec dur = (clock::now() - start);
        fprintf(stderr, "%s%.2fs\n", msg.data(), static_cast<double>(dur.count()));
    }

    std::chrono::time_point<clock> start;
    std::string msg;
};


void print_help() {
#ifndef M64
    fprintf(stderr,
"pfbwt-f. use the prefix-free parsing algorithm to build a BWT for genomic data.\n\
\n"
"usage\n\
    ./pfbwt-f [options] <fasta file>\n\
\n");
#else
    fprintf(stderr,
"pfbwt-f64. use the prefix-free parsing algorithm to build a BWT for BIG genomic data.\n\
\n"
"usage\n\
    ./pfbwt-f64 [options] <fasta file>\n\
\n");
#endif

    fprintf(stderr,
"results\n\
    BWT of input saved to <fasta file>.bwt. Header lines are excluded.\n\
\n\
options\n\
    -s              Build full suffix array and output to <fasta file>.sa\n\
    \n\
    -r              Build run-length sampled suffix arrray and output run-starts to <fasta file>.ssa and run-ends to <fasta file>.esa\n\
    \n\
    -w <int>        window-size for parsing [default: 10] \n\
    \n\
    -p <int>        modulo for parsing [default: 100]\n\
    \n\
    -m              build BWT on external memory\n\
    \n\
    --parse-only    only produce parse (dict, occ, ilist, last, bwlast files), do not build BWT\n\
    \n\
    -h              print this help message\n\
");
}

Args parse_args(int argc, char** argv) {
    Args args;
    int c;

    fputs("==== Command line:", stderr);
    for(int i=0;i<argc;i++)
        fprintf(stderr, " %s",argv[i]);
    fputs("\n", stderr);

    std::string sarg;

    static struct option lopts[] = {
        {"parse-only", no_argument, &args.parse_only, 1}
    };

    while ((c = getopt_long( argc, argv, "w:hsrfm", lopts, NULL) ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 's':
                args.sa = true; break;
            case 'r':
                args.rssa = true; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'm':
                args.mmap = true; break;
            case 'p':
                args.p = atoi(optarg); break;
            case 'h':
                print_help(); exit(0);
            case '?':
                fprintf(stderr, "Unknown option. Use -h for help.\n");
                exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
        args.in_fname.assign( argv[optind] );
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        print_help();
        exit(1);
    }
    return args;
}

template<typename T>
void vec_to_file(const std::vector<T>& vec, std::string fname) {
    FILE* fp = fopen(fname.data(), "wb");
    if (fwrite(vec.data(), sizeof(T), vec.size(), fp) != vec.size() ) {
        die("could not write file");
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

/* saves dict, occs, ilist, bwlast (and bwsai) to disk */
int run_parser(Args args) {
    // build the dictionary and populate .last, .sai and .parse_old
    using parse_t = pfbwtf::Parser<WangHash>;
    parse_t p(args.w, args.p);
    fprintf(stderr, "starting...\n");
    {
        Timer t("TASK\tParsing\t");
        p.parse_fasta(args.in_fname.data(), (args.sa || args.rssa));
    }
    {
        Timer t("TASK\tsorting dict, calculating occs, dumping to file\t");
        FILE* dict_fp = open_aux_file(args.in_fname.data(), EXTDICT, "wb");
        FILE* occ_fp  = open_aux_file(args.in_fname.data(), EXTOCC, "wb");
        p.update_dict([&](const char* phrase, parse_t::UIntType freq) {
                if (fwrite(phrase, 1, strlen(phrase), dict_fp) != strlen(phrase))
                die("Error writing to DICT file\n");
                if (fputc(EndOfWord, dict_fp) == EOF)
                die("Error writing EndOfWord to DICT file");
                if (fwrite(&freq, sizeof(freq), 1, occ_fp) != 1)
                die("Error writing to OCC file\n");
                }
        );
        if (fputc(EndOfDict, dict_fp) == EOF) die("Error writing EndOfDict to DICT file");
        if (fclose(dict_fp)) die("Error closing DICT file");
        else fprintf(stderr, "DICT written to %s.%s\n", args.in_fname.data(), EXTDICT);
        if (fclose(occ_fp)) die("Error closing OCC file");
        else fprintf(stderr, "OCC written to %s.%s\n", args.in_fname.data(), EXTOCC);
    }
    {
        Timer t("TASK\tranking and bwt-ing parse and processing last-chars\t");
        p.bwt_of_parse(
                [&](const std::vector<char>& bwlast, const std::vector<parse_t::UIntType>& ilist, const std::vector<parse_t::UIntType>& bwsai) {
                    vec_to_file<char>(bwlast, args.in_fname + "." + EXTBWLST);
                    vec_to_file<parse_t::UIntType>(ilist, args.in_fname + "." + EXTILIST);
                    if (args.sa || args.rssa) vec_to_file<parse_t::UIntType>(bwsai, args.in_fname + "." + EXTBWSAI);
                },
                (args.sa || args.rssa));
    }
    {
        Timer t("TASK\tdumping files needed by pfbwt\t");
        const auto& parse_ranks = p.get_parse_ranks();
        vec_to_file(parse_ranks, p.get_parse_size(), args.in_fname + "." + EXTPARSE);
    }
    return 0;
}

template<template<typename, typename...> typename R,
         template<typename, typename...> typename W
         >
void run_pfbwt(const Args args) {
    FILE* bwt_fp = open_aux_file(args.in_fname.data(),"bwt","wb");
    using pfbwt_t = pfbwtf::PrefixFreeBWT<R,W>;
    pfbwt_t p(args.in_fname, args.w, args.sa, args.rssa);
    auto bwt_fn = [bwt_fp, args](const char c) {
        fputc(c, bwt_fp);
    };
    if (args.sa) {
        FILE* sa_fp = open_aux_file(args.in_fname.data(), EXTSA, "wb");
        auto sa_fn = [&](const pfbwtf::sa_fn_arg a) {
            fwrite(&a.sa, sizeof(a.sa), 1, sa_fp);
        };
        {
            Timer t("TASK\tgenerating final BWT w/ full SA\t");
            p.generate_bwt_lcp(bwt_fn, sa_fn);
        }
        fclose(sa_fp);
    } else if (args.rssa) {
        FILE* s_fp = open_aux_file(args.in_fname.data(), "ssa", "wb");
        FILE* e_fp = open_aux_file(args.in_fname.data(), "esa", "wb");
        auto sa_fn = [&](const pfbwtf::sa_fn_arg a) {
            switch(a.run_t) {
                case pfbwtf::RunType::START:
                    fwrite(&a.pos, sizeof(a.pos), 1, s_fp);
                    fwrite(&a.sa, sizeof(a.sa), 1, s_fp);
                    break;
                case pfbwtf::RunType::END:
                    fwrite(&a.pos, sizeof(a.pos), 1, e_fp);
                    fwrite(&a.sa, sizeof(a.sa), 1, e_fp);
                    break;
                default:
                    die("error: invalid RunType");
            }
        };
        {
            Timer t("TASK\tgenerating final BWT w/ run-length sampled SA\t");
            p.generate_bwt_lcp(bwt_fn, sa_fn);
        }
        fclose(s_fp);
        fclose(e_fp);
    }
    else { // default case: just output bwt
        {
            Timer t("TASK\tgenerating final BWT w/o SA\t");
            p.generate_bwt_lcp(bwt_fn, [](const pfbwtf::sa_fn_arg a){(void) a;});
        }
    }
    fclose(bwt_fp);
}

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    run_parser(args); // scan file and save relevant info to disk
    if (args.parse_only) return 0;
    if (args.mmap) {
        fprintf(stderr, "workspace will be contained on disk (mmap)\n");
        run_pfbwt<MMapFileSource, MMapFileSink>(args);
    } else {
        fprintf(stderr, "workspace will be contained in memory\n");
        run_pfbwt<VecFileSource, VecFileSinkPrivate>(args);
    }
}
