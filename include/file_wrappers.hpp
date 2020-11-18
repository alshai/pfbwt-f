#ifndef FILEW_HPP
#define FILEW_HPP

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "mio.hpp"

void die_(const char* string) {
    fprintf(stderr, "%s\n", string);
    exit(1);
}

size_t get_file_size_(const char* path) {
    if (!strcmp(path, "-")) {
        die_("cannot get file size from stdin");
    }
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        fprintf(stderr, "%s: ", path);
        die_("error opening file");
    }
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fclose(fp);
    return size;
}

/* load file using mmap and treat it as an STL container (mio)
 * TODO: create a MAP_PRIVATE version
 */

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
        if (!mm.size()) {
            die_("MMapFile: file empty!");
        }
        if (mm.size() %  sizeof(T) != 0) {
            die_("error: file is not evenly divided into size(T)-sized words");
        }
    }

    /* start fresh from a file, given a size
     * size is number of elements, not number of bytes
     */
    void init_file(std::string path, size_t size) {
        // make sure that file is 'cleared' and set to the appropriate size
        FILE* fp = fopen(path.data(), "wb");
        if (fp) {
            fclose(fp);
        } else {
            fprintf(stderr, "%s: ", path.data());
            die_("error opening file");
        }
        if (truncate(path.data(), size * sizeof(T))) {
            die_("error setting file size");
        }
        // map file
        std::error_code e;
        mm.map(path, 0, size * sizeof(T), e);
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

/* NOTE: I think this is MAP_SHARED by default.
 * TODO: add option for MAP_PRIVATE */
template<typename T>
using MMapFileSink = MMapFile<T, mio::access_mode::write>;

/* loads a file into a read-only vector. neither the data nor the underlying
 * file are allowed to be changed */
template<typename T>
class VecFileSource : private std::vector<T> {

    public:

    VecFileSource() = default;

    VecFileSource(std::string path) {
        size_t size = get_file_size_(path.data());
        size_t nelems = size / sizeof(T);
        FILE* fp = fopen(path.data(), "rb");
        if (fp == NULL) {
            fprintf(stderr, "VecFileSource: error opening %s\n", path.data());
            exit(1);
        }
        this->resize(nelems);
        if (fread(&(this->data())[0], sizeof(T), this->size(), fp) != nelems) {
            fprintf(stderr, "VecFileSource: error reading from %s\n", path.data());
            exit(1);
        }
        fclose(fp);
    }

    typename std::vector<T>::const_reference operator[](typename std::vector<T>::size_type i) const {
        return std::vector<T>::operator[](i);
    }

    typename std::vector<T>::size_type size() const {
        return std::vector<T>::size();
    }
};

template<typename T>
class VecFileSink : public std::vector<T> {

    public:

    VecFileSink() = default;

    /* load data from whole file */
    VecFileSink(std::string path) : fname(path) {
        size_t size = get_file_size_(path.data());
        size_t nelems = size / sizeof(T);
        FILE* fp = fopen(path.data(), "rb");
        if (fp == NULL) {
            fprintf(stderr, "VecFileSink: error opening %s\n", path.data());
            exit(1);
        }
        this->resize(nelems);
        if (fread(&(*this)[0], sizeof(T), this->size(), fp) != nelems) {
            fprintf(stderr, "VecFileSink: error reading from %s\n", path.data());
            exit(1);
        }
        fclose(fp);
    }

    /* use when initialized with default constructor in order to reserve heap space
     * for an empty vector and store the file name to which the vector will be written
     * upon destruction
     */
    void init_file(std::string path, size_t s) {
        fname = path;
        this->resize(s);
    }

    protected:

    std::string fname;
};

/* use this when you want allow the underlying file to reflect
 * any changes made to its data */
template<typename T>
class VecFileSinkShared : public VecFileSink<T> {
    public:

    VecFileSinkShared() = default;
    VecFileSinkShared(std::string s) : VecFileSink<T>(s) {};

    /* writes data to file before exiting */
    ~VecFileSinkShared() {
        if (this->fname == "") {
            fprintf(stderr, "no file specified. Data will not be saved to disk\n");
        } else {
            FILE* fp = fopen(this->fname.data(), "wb");
            if (fwrite(&(this->data())[0], sizeof(T), this->size(), fp) != this->size()) {
                die_("unable to write data");
            }
            fclose(fp);
        }
    }
};

template<typename T>
using VecFileSinkPrivate = VecFileSink<T>;
#endif
