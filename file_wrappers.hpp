#ifndef FILEW_HPP
#define FILEW_HPP

#include <vector>
#include <cstdio>
#include <string>
#include "mio.hpp"
extern "C" {
#include "utils.h"
}

/* load file using mmap and treat it as an STL container (mio) */
template<typename T, mio::access_mode AccessMode>
class MMapFile {

    public:

    using value_type = T;
    using size_type = size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator = value_type*;
    using const_iterator = const value_type*;


    MMapFile() = default;

    MMapFile(std::string path) :
        mm(path)
    {
        if (mm.size() %  sizeof(T) != 0) {
            fprintf(stderr, "error: file is not evenly divided into size(T)-sized words\n");
            exit(1);
        }
    }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type>
    reference operator[](const size_type i) noexcept { 
        return reinterpret_cast<T*>(mm.data())[i]; 
    }

    const_reference operator[](const size_type i) const noexcept { 
        return reinterpret_cast<const T*>(mm.data())[i]; 
    }

    size_t size() const noexcept {
        return mm.size() / sizeof(T);
    }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type>
    T* data() noexcept { 
        return reinterpret_cast<T*>(mm.data()); 
    }

    const T* data() const noexcept { 
        return reinterpret_cast<const T*>(mm.data());
    }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type> 
    iterator begin() noexcept { return data(); }
    const_iterator begin() const noexcept { return data(); }
    const_iterator cbegin() const noexcept { return data(); }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type> 
    iterator end() noexcept { return data() + size(); }
    const_iterator end() const noexcept { return data() + size(); }
    const_iterator cend() const noexcept { return data() + size(); }

    private:
    mio::basic_mmap<AccessMode, unsigned char> mm;
    T* data_; // this is derived from mm
};

template<typename T>
using MMapFileSource = MMapFile<T, mio::access_mode::read>;

template<typename T>
using MMapFileSink = MMapFile<T, mio::access_mode::write>;

/* load a file into std::vector using default allocation behaviour (probably malloc) */
template<typename T>
class VecFile : public std::vector<T> {

    public:

    VecFile() {}

    VecFile(std::string path) {
        size_t size = get_file_size(path.data());
        size_t nelems = size / sizeof(T);
        FILE* fp = fopen(path.data(), "rb");
        this->resize(nelems);
        if (fread(&(*this)[0], sizeof(T), this->size(), fp) != nelems) {
            exit(1);
        }
        fclose(fp);
    }

};

#endif
