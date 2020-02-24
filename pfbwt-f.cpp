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
    std::string output;
    std::string stdout_ext;
    size_t w = 10;
    size_t p = 100;
    int sa = 0;
    int rssa = 0;
    int da = 0;
    int ma = 0;
    int mmap = 0;
    int parse_only = 0;
    int trim_non_acgt = 0;
    int non_acgt_to_a = 0;
    int pfbwt_only = 0;
    int verbose = false;
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
    -s                  Output full suffix array to <fasta file>.sa\n\
    \n\
    -r                  Output run-length sampled suffix arrray to \n\
                        <fasta file>.ssa (run-starts) and \n\
                        <fasta file>.esa (run-ends)\n\
    \n\
    -w <int>            window-size for parsing [default: 10] \n\
    \n\
    -p <int>            modulo for parsing [default: 100]\n\
    \n\
    -m                  build BWT on external memory\n\
    \n\
    --parse-only        only produce parse (dict, occ, ilist, last, bwlast)\n\
                        do not build final BWT\n\
    \n\
    -c/--stdout <ext>   output file ending <ext> will be stdout instead.\n\
                        (example: '-c bwt' would output <fasta file>.bwt to stdout)\n\
                        options: bwt, sa\n\
    -h                  print this help message\n\
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
        {"parse-only", no_argument, &args.parse_only, 1},
        {"pfbwt-only", no_argument, &args.pfbwt_only, 1},
        {"trim-non-acgt", no_argument, &args.trim_non_acgt, 1},
        {"non-acgt-to-a", no_argument, &args.non_acgt_to_a, 1},
        {"stdout", required_argument, NULL, 'c'},
        {"verbose", no_argument, &args.verbose, 1},
        {"sa", no_argument, NULL, 's'},
        {"rssa", no_argument, NULL, 'r'},
        {"ma", no_argument, &args.ma, 1},
        {"da", no_argument, NULL, 'd'},
        {"mmap", no_argument, NULL, 'm'},
        {"output", required_argument, NULL, 'o'},
        {"window-size", required_argument, NULL, 'w'},
        {"mod-val", required_argument, NULL, 'p'}
    };

    while ((c = getopt_long( argc, argv, "w:p:o:hsrfm", lopts, NULL) ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 's':
                args.sa = 1; break;
            case 'r':
                args.rssa = 1; break;
            case 'd':
                args.da = 1; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'm':
                args.mmap = 1; break;
            case 'p':
                args.p = atoi(optarg); break;
            case 'h':
                print_help(); exit(0);
            case 'o':
                args.output.assign(optarg); break;
            case 'c':
                args.stdout_ext = optarg; break;
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
        fprintf(stderr, "reading from stdin. Parsing might be a bit slow.\n");
        args.in_fname.assign("-");
    }

    if (args.non_acgt_to_a && args.trim_non_acgt) {
        die("cannot have both --non-acgt-to-a and --trim-non-acgt options enabled at same time");
    }
    if (args.rssa && args.sa) {
        die("cannot have both --sa and --rssa options enabled at same time");
    }

    if (args.in_fname == "-" && args.output == "") {
        die("if reading from stdin, need a prefix for output files (-o, --output)");
    }
    if (args.in_fname != "-"  && args.output == "") {
        args.output = args.in_fname;
    }

    /* TODO: remove this */
    if (args.ma) die("marker array support coming soon!");

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

pfbwtf::ParserParams args_to_parser_params(Args args) {
    pfbwtf::ParserParams p;
    p.fname = args.in_fname;
    p.w = args.w;
    p.p = args.p;
    p.get_sai = args.sa || args.rssa;
    p.get_da = args.da;
    p.verbose = args.verbose;
    p.trim_non_acgt = args.trim_non_acgt;
    p.non_acgt_to_a = args.non_acgt_to_a;
    return p;
}


pfbwtf::PrefixFreeBWTParams args_to_pfbwt_params(Args args) {
    pfbwtf::PrefixFreeBWTParams p;
    p.prefix = args.in_fname;
    p.w = args.w;
    p.sa = args.sa;
    p.rssa  = args.rssa;
    p.ma = args.ma;
    p.da = args.da;
    p.verb = args.verbose;
    return p;
}

/* saves dict, occs, ilist, bwlast (and bwsai) to disk */
int run_parser(Args args) {
    // build the dictionary and populate .last, .sai and .parse_old
    using parse_t = pfbwtf::Parser<WangHash>;
    pfbwtf::ParserParams params(args_to_parser_params(args));
    parse_t p(params);
    fprintf(stderr, "starting...\n");
    {
        Timer t("TASK\tParsing\t");
        p.parse_fasta();
    }
    {
        Timer t("TASK\tsorting dict, calculating occs, dumping to file\t");
        FILE* dict_fp = open_aux_file(args.output.data(), EXTDICT, "wb");
        FILE* occ_fp  = open_aux_file(args.output.data(), EXTOCC, "wb");
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
        else fprintf(stderr, "DICT written to %s.%s\n", args.output.data(), EXTDICT);
        if (fclose(occ_fp)) die("Error closing OCC file");
        else fprintf(stderr, "OCC written to %s.%s\n", args.output.data(), EXTOCC);
    }
    {
        Timer t("TASK\tranking and bwt-ing parse and processing last-chars\t");
        p.bwt_of_parse(
                [&](const std::vector<char>& bwlast, const std::vector<parse_t::UIntType>& ilist, const std::vector<parse_t::UIntType>& bwsai) {
                    vec_to_file<char>(bwlast, args.output + "." + EXTBWLST);
                    vec_to_file<parse_t::UIntType>(ilist, args.output + "." + EXTILIST);
                    if (args.sa || args.rssa) vec_to_file<parse_t::UIntType>(bwsai, args.output + "." + EXTBWSAI);
                    // TODO: should we do DA stuff here too?
                });
    }
    {
        Timer t("TASK\tdumping files needed by pfbwt\t");
        const auto& parse_ranks = p.get_parse_ranks();
        vec_to_file(parse_ranks, p.get_parse_size(), args.output + "." + EXTPARSE);
    }
    if (args.da) {
        const auto& doc_starts = p.get_doc_starts();
        vec_to_file(doc_starts, args.output + ".dstarts");
        const auto& doc_names = p.get_doc_names();
        vec_to_file(doc_names, args.output + ".dnames");
    }
    // TODO: dump ntab to file if applicable.
    if (args.trim_non_acgt) {
        vec_to_file(p.get_ntab(), args.output + ".ntab");
    }
    return 0;
}

std::FILE* init_file_pointer(const Args& args, std::string ext) {
    if (args.stdout_ext == ext) {
        return stdout;
    } else {
        return open_aux_file(args.output.data(), ext.data(), "wb");
    }
}

template<template<typename, typename...> typename R,
         template<typename, typename...> typename W
         >
void run_pfbwt(const Args args) {
    pfbwtf::PrefixFreeBWTParams pfbwt_args(args_to_pfbwt_params(args));
    std::FILE* bwt_fp = init_file_pointer(args, "bwt");
    using pfbwt_t = pfbwtf::PrefixFreeBWT<R,W>;
    pfbwt_t p(pfbwt_args);
    char pc = 0;
    size_t n(0), r(0);
    auto bwt_fn = [&](const char c) {
        fputc(c, bwt_fp);
        if (pc != c) ++r;
        pc = c;
        ++n;
    };
    if (args.sa) {
        std::FILE* sa_fp = init_file_pointer(args, "sa");
        // std::FILE* sa_fp = open_aux_file(args.output.data(), EXTSA, "wb");
        auto sa_fn = [&](const pfbwtf::sa_fn_arg a) {
            fwrite(&a.sa, sizeof(a.sa), 1, sa_fp);
        };
        {
            Timer t("TASK\tgenerating final BWT w/ full SA\t");
            p.generate_bwt_lcp(bwt_fn, sa_fn, [](...){});
        }
        fclose(sa_fp);
    } else if (args.rssa) {
        FILE* s_fp = open_aux_file(args.output.data(), "ssa", "wb");
        FILE* e_fp = open_aux_file(args.output.data(), "esa", "wb");
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
            p.generate_bwt_lcp(bwt_fn, sa_fn, [](...){});
        }
        fclose(s_fp);
        fclose(e_fp);
    }
    else { // default case: just output bwt
        {
            Timer t("TASK\tgenerating final BWT w/o SA\t");
            p.generate_bwt_lcp(bwt_fn, [](const pfbwtf::sa_fn_arg a){(void) a;}, [](...){});
        }
    }
    fprintf(stderr, "n: %lu\n", n);
    fprintf(stderr, "r: %lu\n", r);
    fprintf(stderr, "n/r: %.3f\n", static_cast<double>(n) / r);
    fclose(bwt_fp);
}

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    if (!args.pfbwt_only)
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
