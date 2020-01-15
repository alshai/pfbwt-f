#include <cstdio>
#include <chrono>
#include <getopt.h>
#include "hash.hpp"
#include "parse-f.hpp"
extern "C" {
#include "utils.h"
}

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
    fprintf(stderr, "Usage:\n\t./parse [options] <fasta file>\n");
    fprintf(stderr, "Options:\n-w\twindow size\n-p\tstop word modulus\n-s\tprint SA info\n");
}

pfbwtf::ParseArgs parse_args(int argc, char** argv) {
    pfbwtf::ParseArgs args;
    int c;
    extern char *optarg;
    extern int optind;

    fputs("==== Command line:", stderr);
    for(int i=0;i<argc;i++)
        fprintf(stderr, " %s",argv[i]);
    fputs("\n", stderr);

    std::string sarg;
    while ((c = getopt( argc, argv, "p:w:hvsf") ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 'w':
                sarg.assign( optarg );
                args.w = stoi( sarg ); break;
            case 'p':
                sarg.assign( optarg );
                args.p = stoi( sarg ); break;
            case 's':
                args.sai = true; break;
            case 'h':
                print_help(); exit(1);
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
    // check algorithm parameters
    if ((args.w < 4) || (args.w > 32)) {
        fprintf(stderr, "Windows size must be between 4 and 31 (inclusive)\n");
        exit(1);
    }
    if (args.p<4) {
        fprintf(stderr, "Modulus must be at least 4\n");
        exit(1);
    }
    return args;
}


int main(int argc, char** argv) {
    pfbwtf::ParseArgs args(parse_args(argc, argv));
    // build the dictionary and populate .last, .sai and .parse_old
    pfbwtf::Parser<uint32_t, WangHash> p(args.w, args.p);
    {
        Timer t("TASK\tParsing\t");
        p.parse_fasta(args.in_fname.data(), args.sai);
    }
    {
        Timer t("TASK\tsorting dict, calculating occs\t");
        p.update_dict(args.in_fname.data());
    }
    {
        Timer t("TASK\tRanking parse\t");
        p.generate_parse_ranks();
    }
    {
        Timer t("TASK\tbwt-ing parse and processing last-chars\t");
        p.bwt_of_parse(args.in_fname.data(), false);
    }
    {
        Timer t("TASK\tdumping files needed by pfbwt\t");
        p.dump_parse(args.in_fname.data());
        p.dump_last(args.in_fname.data());
        p.dump_bwlast(args.in_fname.data());
        p.dump_ilist(args.in_fname.data());
    }
    return 0;
}
