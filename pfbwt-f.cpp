#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <getopt.h>
#include "pfbwt-f.hpp"
#include "file_wrappers.hpp"
extern "C" {
#include "utils.h"
}

struct PrefixFreeBWTArgs {
    std::string prefix;
    int w = 10;
    bool sa = false;
    bool rssa = false;
    bool mmap = false;
};

void print_help() {
    return;
}

PrefixFreeBWTArgs parse_args(int argc, char** argv) {
    PrefixFreeBWTArgs args;
    int c;

    fputs("==== Command line:", stderr);
    for(int i=0;i<argc;i++)
        fprintf(stderr, " %s",argv[i]);
    fputs("\n", stderr);

    std::string sarg;
    while ((c = getopt( argc, argv, "w:hsrfm") ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 's':
                args.sa = true; break;
            case 'r':
                args.rssa = true; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'h':
                print_help(); exit(1);
            case 'm':
                args.mmap = true; break;
            case '?':
                fprintf(stderr, "Unknown option. Use -h for help.\n");
                exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
        args.prefix.assign( argv[optind] );
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        print_help();
        exit(1);
    }
    return args;
}

template<typename IntType,
         template<typename, typename...> typename R,
         template<typename, typename...> typename W
         >
void run_pfbwt(const PrefixFreeBWTArgs args) {
    FILE* bwt_fp = open_aux_file(args.prefix.data(),"bwt","wb");
    pfbwtf::PrefixFreeBWT<IntType, R, W> p(args.prefix, args.w, args.sa, args.rssa);
    auto bwt_fn = [bwt_fp, args](const char c) {
        fputc(c, bwt_fp);
    };
    if (args.sa) {
        FILE* sa_fp = open_aux_file(args.prefix.data(), EXTSA, "wb");
        auto sa_fn = [&](const pfbwtf::sa_fn_arg<IntType> a) {
            fwrite(&a.sa, sizeof(IntType), 1, sa_fp);
        };
        p.generate_bwt_lcp(bwt_fn, sa_fn);
        fclose(sa_fp);
    } else if (args.rssa) {
        FILE* s_fp = open_aux_file(args.prefix.data(), "ssa", "wb");
        FILE* e_fp = open_aux_file(args.prefix.data(), "esa", "wb");
        auto sa_fn = [&](const pfbwtf::sa_fn_arg<IntType> a) {
            switch(a.run_t) {
                case pfbwtf::RunType::START:
                    fwrite(&a.pos, sizeof(IntType), 1, s_fp);
                    fwrite(&a.sa, sizeof(IntType), 1, s_fp);
                    break;
                case pfbwtf::RunType::END:
                    fwrite(&a.pos, sizeof(IntType), 1, e_fp);
                    fwrite(&a.sa, sizeof(IntType), 1, e_fp);
                    break;
                default:
                    die("error: invalid RunType");
            }
        };
        p.generate_bwt_lcp(bwt_fn, sa_fn);
        fclose(s_fp);
        fclose(e_fp);
    }
    else { // default case: just output bwt
        p.generate_bwt_lcp(bwt_fn, [](const pfbwtf::sa_fn_arg<IntType> a){(void) a;});
    }
    fclose(bwt_fp);
}

int main(int argc, char** argv) {
    PrefixFreeBWTArgs args(parse_args(argc, argv));
    if (args.mmap) {
        fprintf(stderr, "workspace will be contained on disk (mmap)\n");
        run_pfbwt<uint32_t, MMapFileSource, MMapFileSink>(args);
    } else {
        fprintf(stderr, "workspace will be contained in memory\n");
        run_pfbwt<uint32_t, VecFileSource, VecFileSinkPrivate>(args);
    }
}
