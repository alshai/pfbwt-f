#include <cstdlib>
#include <cinttypes>
#include <cstdio>
#include <string>
#include <iostream>
#include <getopt.h>
#include "pfbwt_io.hpp"
#include "pfparser.hpp"

struct Args {
    std::vector<std::string> prefixes;
    std::string output = "out";
    int w = 10;
    int p = 100; 
    int store_docs = 0;
    // int force = 0;
    // int remove_old_files = 0;
};

Args parse_args(int argc, char** argv) {
    Args args;
    int c;
    static struct option lopts[] = {
        {"docs", no_argument, &args.store_docs, 1},
        {"window-size", required_argument, NULL, 'w'},
        {"mod-val", required_argument, NULL, 'p'},
        {"output", required_argument, NULL, 'o'}
        // {"remove-old-files", no_argument, NULL, 'r'},
        // {"force", no_argument, NULL, 'f'}
    };

    while ((c = getopt_long( argc, argv, "dw:p:o:", lopts, NULL) ) != -1) {
        switch(c) {
            case 'd':
                args.store_docs = 1; break;
            // case 'f':
            //     args.force = 1; break;
            // case 'r':
            //     args.remove_old_files = 1; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'p':
                args.p = atoi(optarg); break;
            case 'o':
                args.output = optarg; break;
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

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    pfbwtf::PfParserParams params;
    params.store_docs = args.docs;
    params.w = args.w;
    params.p = args.p;
    pfbwtf::PfParser<> parser;
    for (auto prefix: args.prefixes) {
        std::cerr << "adding: " << prefix << std::endl;
        parser += pfbwtf::load_parser(prefix, params);
        std::string dict_fname = prefix + ".dict";
        std::string occ_fname = prefix + ".occ";
        std::string parse_ranks_fname = prefix + ".parse";
        // TODO: figures out why this is deleting write-protected files
        // if (args.remove_old_files) {
        //     if (std::remove(dict_fname.data())) fprintf(stderr, "warning: could not remove %s\n", dict_fname.data());
        //     if (std::remove(occ_fname.data())) fprintf(stderr, "warning: could not remove %s\n", occ_fname.data());
        //     if (std::remove(parse_ranks_fname.data())) fprintf(stderr, "warning: could not remove %s\n", parse_ranks_fname.data());
        // }
    }
    parser.finalize();
    pfbwtf::save_parser(parser, args.output);
    return 0;
}
