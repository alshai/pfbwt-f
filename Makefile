# compilation flags
CXX_FLAGS=-std=c++11 -Ofast -Wall -Wextra -march=native -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
CXX=g++

# main executables
EXECS=pfbwt-f64 pfbwt-f

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

vcf_to_bwt: pfbwt-f64 merge_pfp vcf_scan merge_marker_indexes marker_index_to_array


gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

pfbwt-f: pfbwt-f.cpp utils.o gsa/gsacak.o pfbwt.hpp pfparser.hpp file_wrappers.hpp
	$(CXX) $(CXX_FLAGS) -I./ -I./sdsl-lite/include -o $@ pfbwt-f.cpp utils.o gsa/gsacak.o -lz

pfbwt-f64: pfbwt-f.cpp utils.o gsa/gsacak64.o pfbwt.hpp pfparser.hpp file_wrappers.hpp pfbwt_io.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -I./ -I./sdsl-lite/include -o $@ pfbwt-f.cpp utils.o gsa/gsacak64.o -lz

simplebwt: simplebwt.o gsa/gsacak.o
	$(CC) $(CFLAGS) -o $@ $< gsa/gsacak.o

simplebwt64: simplebwt.o gsa/gsacak64.o
	$(CC) $(CFLAGS) -DM64 -o $@ $< gsa/gsacak64.o

dump_intfile: scripts/dump_intfile.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $<

test_parser: pfparser.hpp pfbwt_io.hpp tests/test_parser.cpp utils.o
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ tests/test_parser.cpp utils.o -lz

merge_pfp: merge_pfp.cpp pfparser.hpp pfbwt_io.hpp utils.o
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ merge_pfp.cpp gsa/gsacak64.o utils.o -lz -lpthread

vcf_scan: marker_array/vcf_scan.cpp marker_array/vcf_scanner.hpp marker_array/marker_index.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/vcf_scan.cpp -lhts -lsdsl

generate_marker_array: marker_array/generate_marker_array.cpp marker_array/marker_index.hpp marker_array/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/generate_marker_array.cpp -lsdsl

merge_marker_indexes: marker_array/merge_marker_indexes.cpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/merge_marker_indexes.cpp

load_marker_index: marker_array/load_marker_index.cpp marker_array/marker_index.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/load_marker_index.cpp utils.o -lsdsl

marker_index_to_array: marker_array/marker_index_to_array.cpp marker_array/marker_index.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/marker_index_to_array.cpp utils.o -lsdsl

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

test: quicktest

quicktest:
	python tests/quick_test.py --sa

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
