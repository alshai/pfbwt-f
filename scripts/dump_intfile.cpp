#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <getopt.h>
#include <type_traits>

template<typename T>
void dump_intfile(std::FILE* fp) {
    T x;
    while (fread(&x, sizeof(T), 1, fp) == 1) {
        fprintf(stdout, "%lu\n", static_cast<uint64_t>(x)); // TODO: customize based on inttype
    }
}

enum Bytes {THIRTYTWO, SIXTYFOUR}; // , STHIRTYTWO, SSIXTYFOUR};

int main(int argc, char** argv) {
    int c;
    Bytes b = THIRTYTWO;
    while ((c = getopt( argc, argv, "bl") ) != -1) {
        switch(c) {
            case 'b':
               b = THIRTYTWO; break;
            case 'l':
                b = SIXTYFOUR; break;
            case '?':
                fprintf(stderr, "Unknown option. Use -h for help.\n");
                exit(1);
        }
    }
    // the only input parameter is the file name
    std::FILE* fp = NULL;
    if (argc == optind+1) {
        fp = std::fopen(argv[optind], "rb");
    }
    else {
        fp = std::freopen(NULL, "rb", stdin);
    }

    switch (b) {
        case THIRTYTWO:
            dump_intfile<uint32_t>(fp); break;
        case SIXTYFOUR:
            dump_intfile<uint64_t>(fp); break;
    }
    fclose(fp);
}