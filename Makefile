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

vcf_to_bwt: pfbwt-f64 merge_pfp vcf_scan merge_mps mps_to_ma


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

vcf_scan: marker_array/vcf_scan.cpp marker_array/vcf_scanner.hpp marker_array/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/vcf_scan.cpp -lhts -I./ -I./sdsl-lite/include

generate_marker_array: marker_array/generate_marker_array.cpp marker_array/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/generate_marker_array.cpp -I./sdsl-lite/include

merge_mps: marker_array/merge_mps.cpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/merge_mps.cpp

load_rlarr: marker_array/load_rlarr.cpp marker_array/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/load_rlarr.cpp utils.o -I./sdsl-lite/include

dump_markers: marker_array/dump_markers.cpp marker_array/marker.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/dump_markers.cpp

mps_to_ma: marker_array/mps_to_ma.cpp marker_array/marker_array.hpp marker_array/rle_window_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ marker_array/mps_to_ma.cpp utils.o -I./ -I./sdsl-lite/include

test_load: tests/test_load.cpp marker_array/marker_array.hpp marker_array/rle_window_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ tests/test_load.cpp -I./ -I./sdsl-lite/include

test_markers: tests/test_markers.cpp marker_array/marker_array.hpp marker_array/rle_window_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ tests/test_markers.cpp -I./ -I./sdsl-lite/include

marker_array/rle_window_array.hpp: sdsl_bv_wrappers.hpp file_wrappers.hpp

marker_array/marker_array.hpp: marker_array/marker_array.hpp marker_array/rle_window_array.hpp file_wrappers.hpp marker_array/marker.hpp

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

test: quicktest

quicktest:
	python tests/quick_test.py --sa

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
