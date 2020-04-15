#include <cstdlib>
#include <cinttypes>
#include <cstdio>
#include <string>
#include <iostream>
#include <getopt.h>
#include <thread>
#include "pfbwt_io.hpp"
#include "pfparser.hpp"

struct Args {
    std::vector<std::string> prefixes;
    std::string output = "out";
    size_t nthreads = 1;
    int w = 10;
    int p = 100;
    int store_docs = 0;
};

Args parse_args(int argc, char** argv) {
    Args args;
    int c;
    static struct option lopts[] = {
        {"docs", no_argument, &args.store_docs, 1},
        {"window-size", required_argument, NULL, 'w'},
        {"mod-val", required_argument, NULL, 'p'},
        {"output", required_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'}
    };

    while ((c = getopt_long( argc, argv, "dw:p:o:t:", lopts, NULL) ) != -1) {
        switch(c) {
            case 'd':
                args.store_docs = 1; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'p':
                args.p = atoi(optarg); break;
            case 'o':
                args.output = optarg; break;
            case 't':
                args.nthreads = atoi(optarg); break;
            case '?':
                std::cerr << "Unknown option.\n";
                exit(1);
                break;
            case ':':
                std::cerr << "no argument specified for option\n";
                exit(1);
        }
    }
    for (int i = optind; i < argc; ++i) {
        args.prefixes.push_back(argv[optind++]);
    }
    if (args.output == "") {
        std::cerr << "no output prefix specified. Defaulting to 'out'\n";
    }

    return args;
}

using ParseV = std::vector<pfbwtf::PfParser<>>;

struct MergeArgs {
    ~MergeArgs() { for (auto fp: logs) if (fp != NULL) fclose(fp); }
    void init_logs() {
        logs.clear();
        for (size_t i = 0; i < nthreads; ++i) {
            std::string fname = output + ".pfbwtf.th_" + std::to_string(i) + ".log";
            logs.push_back(fopen(fname.data(), "w"));
            if (logs.back() == NULL) exit(1);
        }
    }
    std::vector<std::string> prefixes;
    std::vector<FILE*> logs;
    std::string output = "out";
    pfbwtf::PfParserParams params;
    size_t nthreads;
    ParseV parsers;
};

void parser_merge_worker(MergeArgs& args, size_t tidx, size_t start_i, size_t end_i) {
    for (size_t i = start_i; i <= end_i; ++i) {
        std::string prefix = args.prefixes[i];
        args.parsers[tidx] += pfbwtf::load_or_generate_parser_w_log(prefix, args.params, args.logs[tidx]);
    }
    args.parsers[tidx].finalize();
}

pfbwtf::PfParser<> parser_merge_from_vec(ParseV& pv) {
    pfbwtf::PfParser<> parser;
    for (const auto& p: pv) {
        parser += p;
    }
    parser.finalize();
    return parser;
}

void merge_pfp(Args args) {
    pfbwtf::PfParserParams params;
    params.store_docs = args.store_docs;
    params.w = args.w;
    params.p = args.p;
    if (args.nthreads > 1) {
        // initialize threads and thread arguments
        MergeArgs margs;
        margs.prefixes = args.prefixes;
        margs.nthreads = args.nthreads;
        margs.params = params;
        margs.output = args.output;
        margs.init_logs();
        margs.parsers.resize(args.nthreads);
        std::vector<std::thread> threads;
        threads.reserve(args.nthreads);
        for (size_t i = 0; i < args.nthreads; ++i) {
            size_t start = i * args.prefixes.size() / args.nthreads;
            size_t end = ((i + 1) * args.prefixes.size() / args.nthreads) - 1;
            end = end > args.prefixes.size() - 1 ? args.prefixes.size() - 1 : end;
            threads.push_back(std::thread(parser_merge_worker, std::ref(margs), i, start, end));
        }
        for (size_t i = 0; i < args.nthreads; ++i)  {
            threads[i].join();
        }
        auto parser = parser_merge_from_vec(margs.parsers);
        pfbwtf::save_parser(parser, args.output);
    } else {
        std::string log_fname = args.output + ".pfbwt.log";
        FILE* fp = fopen(log_fname.data(), "w");
        if (fp == NULL) {fprintf(stderr, "error opening log\n"); exit(1);}
        pfbwtf::PfParser<> parser;
        for (auto prefix: args.prefixes) {
            parser += pfbwtf::load_or_generate_parser_w_log(prefix, params, fp);
        }
        parser.finalize();
        pfbwtf::save_parser(parser, args.output);
    }
}

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    merge_pfp(args);
    return 0;
}
