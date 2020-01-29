# pfbwt-f

Implementation of the Prefix-Free BWT (pfbwt) algorithm desribed in [Boucher, et. al. 2019](https://doi.org/10.1186/s13015-019-0148-5),
optimized specifically for **fasta** files containing genomic data (ie., text
contains only the characters A, C, G, and T. Some Ns are allowed).

Borrows heavily from the [original implementation by Giovanni Manzini](https://gitlab.com/manzai/Big-BWT/)

Some modifications / additions from the original code:

* Compatibility with C++11

* Parsing and BWT-construction steps are now callable from C++11 code.

* Uses the [Wang hash](http://www.burtleburtle.net/bob/hash/integer.html)
instead of a rolling Rabin-Karp hash to select phrases in the parsing step.
(need to be careful when there are a lot of Ns in the input)

* Modifies some data structures during parse step for modest time/memory savings

* Option to use mmap (memory mapping) on workspace data. Allows for BWT to be
constructed **even when the size of the workspace needed for this algorithm exceeds RAM capacity**!

Features that will be added soon.

* ~~Full,~~ Sampled and/or Run-Length Suffix Array.
* Document Array and other SA-related data structures.
* Thread support
* Clearer separation of 32-bit and 64-bit modes

## Installation

To download:

```
git clone --recursive https://github.com/alshai/pfbwt-f
cd pfbwt-f
make
```

## Usage

```
./pfbwt-f -p <mod value> -w <window size> <x.fa>
```

The final BWT will be located in `x.fa.bwt`. This file will contain one extra
character from the input (0x00) - this represents the end-of-string symbol of
the text.

Please use `pfbwt-f64` if your data exceeds 2^32 characters, otherwise results will be incorrect.

## Some features

Output the full Suffix Array to `<x.fa>.sa`:

```
./pfbwt-f -s <x.fa>
```

The output is formatted as consecutive 32-bit integers, no delimiters (64-bit support coming soon!).

Output the run-length Suffix Array to `<x.fa>.ssa` (run-starts) and `<x.fa>.esa` (run-ends):

```
./pfbwt-f -r <x.fa>
```

The output is formatted as consecutive **pairs** of 32-bit integers (position in BWT, SA sample), no delimiters (64-bit support coming soon!).

## Options

```
pfbwt-f. use the prefix-free parsing algorithm to build a BWT for genomic data.

usage
    ./pfbwt-f [options] <fasta file>

    results
        BWT of input saved to <fasta file>.bwt. Header lines are excluded.

    options
        -s              Build full suffix array and output to <fasta file>.sa

        -r              Build run-length sampled suffix arrray and output run-starts to <fasta file>.ssa and run-ends to <fasta file>.esa

        -w <int>        window-size for parsing [default: 10]

        -p <int>        modulo for parsing [default: 100]

        -m              build BWT on external memory

        --parse-only    only produce parse (dict, occ, ilist, last, bwlast files), do not build BWT

        -h              print this help message
```
