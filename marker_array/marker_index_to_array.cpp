#include <string>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include "marker_array/marker_index.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: ./marker_index_to_array <marker index> <sa> [<output prefix>] \n");
    }
    if (argc > 3) {
        write_marker_array(argv[1], argv[2], argv[3]);
    } else {
        write_marker_array(argv[1], argv[2]);
    }
    return 0;
}
