# Copyright (C) 2015 Diego Darriba
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <Diego.Darriba@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

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
CPPFLAGS = -std=c++11
LDFLAGS = -lpll_tree -lpll -lm -lsdsl -ldivsufsort -ldivsufsort64 -lstdc++

OBJS = main.o modified_library_functions.o util.o 
PROG = main

default: all
all : $(PROG)

main : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

%.o: %.c %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) $(PROG)
