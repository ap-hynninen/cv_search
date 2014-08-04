
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

SRC = cv.cpp cv_util.cpp LM.cpp

OBJS = cv.o cv_util.o LM.o

ifeq ($(OS),Linux)
LFLAGS = -std=c++0x
else
LFLAGS = -lstdc++.6
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
	$(CC) -c -O3 -std=c++11 -D$(DEFS) $<

# DO NOT DELETE

cv.o: LM.hpp cv_util.h
LM.o: LM.hpp
