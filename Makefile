
# Detect OS
OS := $(shell uname -s)

# Detect Intel compiler
INTEL_COMPILER := $(shell which icc | wc -l 2> /dev/null)

DEFS = USE_OPENMP

ifeq ($(INTEL_COMPILER),1)
CC = icc
CL = icc
else
CC = gcc
CL = gcc
endif

SRC = cv.cpp cv_util.cpp LM.cpp Coord.cpp

OBJS = cv.o cv_util.o LM.o Coord.o

CFLAGS = -O0 -g
LFLAGS = -O0 -g

ifeq ($(DEFS),USE_OPENMP)
ifeq ($(INTEL_COMPILER),1)
CFLAGS += -openmp
LFLAGS += -openmp
else
CFLAGS += -fopenmp
LFLAGS += -fopenmp
endif
endif

ifeq ($(INTEL_COMPILER),1)
LFLAGS += -std=c++0x
else
LFLAGS += -lstdc++.6
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
