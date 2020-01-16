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


template<typename T>
void vec_to_file(const std::vector<T>& vec, std::string fname) {
    FILE* fp = fopen(fname.data(), "wb");
    if (fwrite(vec.data(), sizeof(T), vec.size(), fp) != vec.size() ) {
        die("could not write file");
    }
}

template<typename T>
void vec_to_file(const std::vector<T>& vec, size_t nelems, std::string fname) {
    FILE* fp = fopen(fname.data(), "wb");
    if (fwrite(vec.data(), sizeof(T), nelems, fp) != nelems ) {
        die("could not write file");
    }
}

int main(int argc, char** argv) {
    pfbwtf::ParseArgs args(parse_args(argc, argv));
    // build the dictionary and populate .last, .sai and .parse_old
    pfbwtf::Parser<uint32_t, WangHash> p(args.w, args.p);
    fprintf(stderr, "starting...\n");
    {
        Timer t("TASK\tParsing\t");
        p.parse_fasta(args.in_fname.data(), args.sai);
    }
    {
        Timer t("TASK\tsorting dict, calculating occs, dumping to file\t");
        FILE* dict_fp = open_aux_file(args.in_fname.data(), EXTDICT, "wb");
        FILE* occ_fp  = open_aux_file(args.in_fname.data(), EXTOCC, "wb");
        p.update_dict([&](const char* phrase, uint32_t freq) {
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
                [&](std::vector<char> bwlast, std::vector<uint32_t> ilist, std::vector<uint32_t> bwsai) {
                    vec_to_file<char>(bwlast, args.in_fname + "." + EXTBWLST);
                    vec_to_file<uint32_t>(ilist, args.in_fname + "." + EXTILIST);
                    if (args.sai) vec_to_file<uint32_t>(bwsai, args.in_fname + "." + EXTBWSAI);
                }, 
                args.sai);
    }
    {
        Timer t("TASK\tdumping files needed by pfbwt\t");
        const auto& parse_ranks = p.get_parse_ranks();
        vec_to_file(parse_ranks, p.get_parse_size(), args.in_fname + "." + EXTPARSE);
    }
    return 0;
}
