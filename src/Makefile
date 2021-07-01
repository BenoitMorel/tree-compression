# Profiling
OS = $(shell uname -s)
PROFILING=-g
CC ?= clang

# Compiler warnings
ifeq ($(CC),clang)
	ADD_WARN=-Weverything -Wno-padded -Wno-missing-noreturn \
	         -Wno-format-nonliteral -Wno-unused-parameter
endif
WARN=-Wall -Wsign-compare $(ADD_WARN)

CFLAGS = -g -O3 -Wall -Wsign-compare $(PROFILING) $(WARN)
CPPFLAGS = -std=c++11 -I./deps/
LDFLAGS = -L./deps/ -lpllmodtree -lpllmodutil -lpllmodalgorithm -lpllmodbinary -lpllmodmsa -lpllmodoptimize -lpll -lm -lsdsl -ldivsufsort -ldivsufsort64 -lstdc++
INCLUDEFLAGS = -Ldeps/

OBJS = main.o modified_library_functions.o util.o compress_functions.o uncompress_functions.o datastructure_compression_functions.o
PROG = main

default: all
all : $(PROG)

main : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

%.o: %.c %.cpp
	$(CC) $(CPPFLAGS) $(INCLUDEFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) $(PROG)
