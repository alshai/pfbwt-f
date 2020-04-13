#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <getopt.h>
#include <chrono>
#include "pfbwt.hpp"
#include "pfparser.hpp"
#include "hash.hpp"
#include "file_wrappers.hpp"
#include "pfbwt_io.hpp"
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
    int mmap = 0;
    int parse_only = 0;
    int trim_non_acgt = 0;
    int non_acgt_to_a = 0;
    int pfbwt_only = 0;
    int verbose = false;
    int print_docs = 0;
    size_t n = 0;
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
        {"print-docs", no_argument, &args.print_docs, 1},
        {"stdout", required_argument, NULL, 'c'},
        {"verbose", no_argument, &args.verbose, 1},
        {"sa", no_argument, NULL, 's'},
        {"rssa", no_argument, NULL, 'r'},
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
    /*
    if (args.rssa && args.sa) {
        die("cannot have both --sa and --rssa options enabled at same time");
    }
    */

    if (args.in_fname == "-" && args.output == "") {
        die("if reading from stdin, need a prefix for output files (-o, --output)");
    }
    if (args.in_fname != "-"  && args.output == "") {
        args.output = args.in_fname;
    }

    return args;
}


pfbwtf::PfParserParams args_to_parser_params(Args args) {
    pfbwtf::PfParserParams p;
    p.w = args.w;
    p.p = args.p;
    p.get_sai = args.sa || args.rssa;
    p.verbose = args.verbose;
    p.trim_non_acgt = args.trim_non_acgt;
    p.non_acgt_to_a = args.non_acgt_to_a;
    p.store_docs = args.print_docs;
    return p;
}


pfbwtf::PrefixFreeBWTParams args_to_pfbwt_params(Args args) {
    pfbwtf::PrefixFreeBWTParams p;
    p.prefix = args.output;
    p.w = args.w;
    p.sa = args.sa;
    p.rssa  = args.rssa;
    p.verb = args.verbose;
    return p;
}

/* saves dict, occs, ilist, bwlast (and bwsai) to disk */
size_t run_parser(Args args) {
    // build the dictionary and populate .last, .sai and .parse_old
    using parse_t = pfbwtf::PfParser<WangHash>;
    size_t n = 0;
    pfbwtf::PfParserParams params(args_to_parser_params(args));
    parse_t p(params);
    fprintf(stderr, "starting...\n");
    { // TODO: add option for more fasta files
        Timer t("TASK\tparsing input\t");
        n = p.add_fasta(args.in_fname);
    }
    {
        Timer t("TASK\tfinalizing parse, writing dict, occs, and ranks\t");
        p.finalize();
        save_parser(p, args.output);
    }
    fprintf(stderr, "XXX\n");
    {
        Timer t("TASK\tranking and bwt-ing parse and processing last-chars\t");
        p.bwt_of_parse(
                [&](const std::vector<char>& bwlast,
                    const std::vector<parse_t::UIntType>& ilist,
                    const std::vector<parse_t::UIntType>& bwsai) {
                    vec_to_file<char>(bwlast, args.output + "." + EXTBWLST);
                    vec_to_file<parse_t::UIntType>(ilist, args.output + "." + EXTILIST);
                    if (args.sa || args.rssa) vec_to_file<parse_t::UIntType>(bwsai, args.output + "." + EXTBWSAI);
                });
    }
    if (args.print_docs) {
        std::FILE* doc_fp = open_aux_file(args.output.data(), "docs", "w");
        const auto& doc_names = p.get_doc_names();
        const auto& doc_starts = p.get_doc_starts();
        for (size_t i = 0; i < doc_starts.size(); ++i) {
            fprintf(doc_fp, "%s %lu\n", doc_names[i].data(), static_cast<uint64_t>(doc_starts[i]));
        }
        fclose(doc_fp);
    }
    // TODO: dump ntab to file if applicable.
    if (args.trim_non_acgt) {
        vec_to_file(p.get_ntab(), args.output + ".ntab");
    }
    std::FILE* n_fp = open_aux_file(args.output.data(), "n", "w");
    fprintf(n_fp, "%lu\n", n);
    fclose(n_fp);
    return n;
}

std::FILE* init_file_pointer_wb(const Args& args, std::string ext) {
    if (args.stdout_ext == ext) {
        return stdout;
    } else {
        return open_aux_file(args.output.data(), ext.data(), "wb");
    }
}

size_t read_single_int_str(const char* fname, const char* ext) {
    size_t n = 0;
    std::FILE* n_fp = open_aux_file(fname, ext, "r");
    char* line = NULL;
    size_t x = 0;
    if (getline(&line, &x, n_fp) != -1) {
        n = std::atol(line);
    } else {
        die("could not read '.n' file");
        return 0;
    }
    free(line);
    fclose(n_fp);
    return n;
}

template<template<typename, typename...> typename R,
         template<typename, typename...> typename W
         >
void run_pfbwt(const Args args) {
    pfbwtf::PrefixFreeBWTParams pfbwt_args(args_to_pfbwt_params(args));
    std::FILE* bwt_fp = init_file_pointer_wb(args, "bwt");
    using pfbwt_t = pfbwtf::PrefixFreeBWT<R,W>;
    pfbwt_t p(pfbwt_args);
    size_t r = 0;
    size_t n = args.n;
    if (!args.n) {
        fprintf(stderr, "reading n from file\n");
        n = read_single_int_str(args.output.data(), "n");
    }
    if (args.sa | args.rssa ) {
        std::FILE* sa_fp = NULL;
        std::FILE* ssa_fp = NULL;
        std::FILE* esa_fp = NULL;
        if (args.sa)
            sa_fp = init_file_pointer_wb(args, "sa");
        if (args.rssa) {
            ssa_fp = open_aux_file(args.output.data(), "ssa", "wb");
            esa_fp = open_aux_file(args.output.data(), "esa", "wb");
        }
        typename pfbwt_t::UIntType psa = 0;
        typename pfbwt_t::UIntType pi = 0, i = 0;
        auto out_fn = [&](const pfbwtf::out_fn_arg a) {
            fwrite(&a.bwtc, sizeof(a.bwtc), 1, bwt_fp);
            if (args.sa) {
                typename pfbwt_t::UIntType x = i ? a.sa : n;
                fwrite(&x, sizeof(x), 1, sa_fp);
            }
            if (a.bwtc != a.pbwtc) { // run_start
                ++r;
                if (args.rssa) {
                    typename pfbwt_t::UIntType x = i ? a.sa : n;
                    fwrite(&i, sizeof(i), 1, ssa_fp);
                    fwrite(&x, sizeof(x), 1, ssa_fp);
                    if (i) {
                        typename pfbwt_t::UIntType y = pi ? psa : n;
                        fwrite(&pi, sizeof(pi), 1, esa_fp);
                        fwrite(&y, sizeof(y), 1, esa_fp);
                    }
                }
            }
            pi = i;
            psa = a.sa;
            i += 1;
        };
        {
            Timer t("TASK\tgenerating final BWT w/ full and/or run-length SA\t");
            p.generate_bwt_lcp(out_fn);
            // write final run
            if (args.rssa) {
                fwrite(&pi, sizeof(pi), 1, esa_fp);
                fwrite(&psa, sizeof(psa), 1, esa_fp);
            }
        }
        if (args.sa) fclose(sa_fp);
        if (args.rssa) {
            fclose(ssa_fp);
            fclose(esa_fp);
        }
    } else { // default case: just output bwt
        {
            auto out_fn = [&](const pfbwtf::out_fn_arg a) {
                if (a.bwtc != a.pbwtc) ++r;
                fwrite(&a.bwtc, sizeof(a.bwtc), 1, bwt_fp);
            };
            Timer t("TASK\tgenerating final BWT w/o SA\t");
            p.generate_bwt_lcp(out_fn);
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
        args.n = run_parser(args); // scan file and save relevant info to disk
    if (args.parse_only) return 0;
    if (args.mmap) {
        fprintf(stderr, "workspace will be contained on disk (mmap)\n");
        run_pfbwt<MMapFileSource, MMapFileSink>(args);
    } else {
        fprintf(stderr, "workspace will be contained in memory\n");
        run_pfbwt<VecFileSource, VecFileSinkPrivate>(args);
    }
}
