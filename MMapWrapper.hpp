#ifndef MMAPW_HPP
#define MMAPW_HPP

#include <vector>
#include <cstdio>
#include <string>
#include "mio.hpp"

// READ-ONLY abstraction of mio::ummap_source that treats underlying file as
// contiguous array of sizeof(T) byte words
// TODO: add iterators, add write-support
template<typename T, mio::access_mode AccessMode>
class MMapWrapper {

    public:

    using value_type = T;
    using size_type = size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator = value_type*;
    using const_iterator = const value_type*;


    MMapWrapper() = default;

    MMapWrapper(std::string path) :
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

    size_t size() const {
        return mm.size() / sizeof(T);
    }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type>
    T* data() { 
        return reinterpret_cast<T*>(mm.data()); 
    }

    const T* data() const { 
        return reinterpret_cast<const T*>(mm.data());
    }

    /* TODO: to implement, write length() function
    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type> 
    iterator begin() noexcept { return data(); }
    const_iterator begin() const noexcept { return data(); }
    const_iterator cbegin() const noexcept { return data(); }

    template<mio::access_mode A = AccessMode, typename = typename std::enable_if<A == mio::access_mode::write>::type> 
    iterator end() noexcept { return data() + length(); }
    const_iterator end() const noexcept { return data() + length(); }
    const_iterator cend() const noexcept { return data() + length(); }
    */

    private:
    mio::basic_mmap<AccessMode, unsigned char> mm;
    T* data_; // this is derived from mm
};

template<typename T>
using MMapSource = MMapWrapper<T, mio::access_mode::read>;

template<typename T>
using MMapSink = MMapWrapper<T, mio::access_mode::write>;

#endif
