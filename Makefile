# compilation flags
CXX_FLAGS=-std=c++11 -Ofast -Wall -Wextra -g -march=native
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
CXX=g++

# main executables
EXECS=pfbwt-f64 pfbwt-f

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

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
	$(CXX) $(CXX_FLAGS) -DM64 -g -o $@ tests/test_parser.cpp utils.o -lz

merge_pfp: merge_pfp.cpp pfparser.hpp pfbwt_io.hpp utils.o
	$(CXX) $(CXX_FLAGS) -DM64 -g -o $@ merge_pfp.cpp utils.o -lz

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

test: quicktest

quicktest:
	python tests/quick_test.py --sa

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
