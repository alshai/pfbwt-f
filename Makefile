# compilation flags
CXX_FLAGS=-std=c++11 -Ofast -Wall -Wextra -march=native -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
CXX=g++
INC=-I./include
SDSL_INC=-I./sdsl-lite/include

# main executables
EXECS=pfbwt-f64 merge_pfp vcf_scan merge_mps mps_to_ma

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

vcf_to_bwt: $(EXECS)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

pfbwt-f: src/pfbwt-f.cpp src/utils.o gsa/gsacak.o include/pfbwt.hpp include/pfparser.hpp include/file_wrappers.hpp
	$(CXX) $(CXX_FLAGS)  -o $@ src/pfbwt-f.cpp src/utils.o gsa/gsacak.o -lz -I./sdsl-lite/include $(INC)

pfbwt-f64: src/pfbwt-f.cpp src/utils.o gsa/gsacak64.o include/pfbwt.hpp include/pfparser.hpp include/file_wrappers.hpp include/pfbwt_io.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/pfbwt-f.cpp src/utils.o gsa/gsacak64.o -lz $(INC) $(SDSL_INC)

dump_intfile: scripts/dump_intfile.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $<

merge_pfp: src/merge_pfp.cpp include/pfparser.hpp include/pfbwt_io.hpp src/utils.o
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/merge_pfp.cpp gsa/gsacak64.o src/utils.o -lz -lpthread $(INC)

vcf_scan: src/vcf_scan.cpp include/vcf_scanner.hpp include/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/vcf_scan.cpp -lhts $(INC) $(SDSL_INC)

merge_mps: src/merge_mps.cpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/merge_mps.cpp $(INC)

load_rlarr: src/load_rlarr.cpp include/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/load_rlarr.cpp src/utils.o $(INC) $(SDSL_INC)

dump_markers: src/dump_markers.cpp include/marker.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/dump_markers.cpp $(INC)

mps_to_ma: src/mps_to_ma.cpp include/marker_array.hpp
	$(CXX) $(CXX_FLAGS) -DM64 -o $@ src/mps_to_ma.cpp src/utils.o $(INC) $(SDSL_INC)

src/utils.o: src/utils.c include/utils.h
	$(CC) $(CFLAGS) -c -o $@ $< $(INC)

test: quicktest

quicktest:
	python tests/quick_test.py --sa

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
