#include <string>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <getopt.h>
#include "marker_array.hpp"
#include "file_wrappers.hpp"

struct Args {
    std::string mai_fname = "";
    std::string sa_fname = "";
    std::string output = "out";
    int mmap = 0;
};

Args parse_args(int argc, char** argv) {
    Args args;
    int c;
    static struct option lopts[] = {
        {"mmap", no_argument, NULL, 'm'},
        {"output", required_argument, NULL, 'o'}
    };
    while ((c = getopt_long( argc, argv, "o:mh", lopts, NULL) ) != -1) {
        switch(c) {
            case 'm':
                args.mmap = 1; break;
            case 'o':
                args.output = std::string(optarg); break;
            case '?':
                fprintf(stderr,  "Unknown option.\n");
                exit(1);
                break;
            case ':':
                fprintf(stderr,  "no argument specified for option\n");
                exit(1);
        }
    }
    args.mai_fname = argv[optind++];
    args.sa_fname = argv[optind++];
    return args;
}

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    if (args.mmap) {
        write_marker_array<MarkerPositions<MMapFileSource>>(args.mai_fname, args.sa_fname, args.output);
    } else {
        write_marker_array<MarkerPositions<VecFileSource>>(args.mai_fname, args.sa_fname, args.output);
    }
    return 0;
}
