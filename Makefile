
# Detect OS
OS := $(shell uname -s)

# Detect Intel compiler
INTEL_COMPILER := $(shell which icc | wc -l 2> /dev/null)

DEFS = DUMMY

ifeq ($(INTEL_COMPILER),1)
CC = icc
CL = icc
else
CC = gcc
CL = gcc
endif

SRC = cv.cpp cv_util.cpp LM.cpp Coord.cpp

OBJS = cv.o cv_util.o LM.o Coord.o

ifeq ($(INTEL_COMPILER),1)
CFLAGS = -O3 -g -openmp
LFLAGS = -O3 -g -std=c++0x -openmp
else
CFLAGS = -O3 -g -fopenmp
LFLAGS = -O3 -g -lstdc++.6 -fopenmp
endif

all: cv

cv : $(OBJS)
	$(CL) $(LFLAGS) -o cv $(OBJS)

clean: 
	rm -f *.o
	rm -f *~
	rm -f cv

depend:
	makedepend $(SRC)

%.o : %.cpp
	$(CC) -c $(CFLAGS) -D$(DEFS) $<

# DO NOT DELETE

cv.o: LM.hpp cv_util.h
LM.o: LM.hpp
Coord.o: Coord.hpp
