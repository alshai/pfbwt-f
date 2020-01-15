#ifndef VECW_HPP
#define VECW_HPP

#include <cstdio>
#include <string>
#include <vector>

// READ-ONLY wrapper of std::vector such that the constructor takes in a file
// name.  this wrapper exists so that we have a unified interface between using
// std::vector (i.e. malloc) and mio::basic_mmap
template<typename T>
class VecWrapper : public std::vector<T> {

    public:

    VecWrapper() {}

    VecWrapper(std::string path) {
        FILE* fp = fopen(path.data(), "rb");
        fseek(fp, 0, SEEK_END);
        size_t size = ftell(fp);
        size_t nelems = size / sizeof(T);
        rewind(fp);
        this->resize(nelems);
        if (fread(&(*this)[0], sizeof(T), this->size(), fp) != nelems) {
            exit(1);
        }
        fclose(fp);
    }

};

#endif
