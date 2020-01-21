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

1) Parse the input

```
./parse-f -w <window size> -p <hash mod value> <x.fa>
```

NOTE: The input file should not contain the characters 0x00, 0x01, and 0x02,
each of which are used by the internal algorithms.

2) build the BWT. Make sure the window size is the same one as used in `parse-f`.

```
./pfbwt-f -w <window size> <x.fa>
```

The final BWT will be located in `x.fa.bwt`. This file will contain one extra
character from the input (0x00) - this represents the end-of-string symbol of
the text.
