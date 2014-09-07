
# Detect OS
OS := $(shell uname -s)

YES := $(shell which make | wc -l 2> /dev/null)

# Detect Intel compiler
INTEL_COMPILER := $(shell which icc | wc -l 2> /dev/null)

DEFS = -DUSE_OPENMP

ifeq ($(INTEL_COMPILER), $(YES))
CC = icc
CL = icc
else
CC = gcc
CL = gcc
endif

CFLAGS = -O3 
LFLAGS = -O3

ifeq ($(DEFS),-DUSE_OPENMP)
ifeq ($(INTEL_COMPILER), $(YES))
CFLAGS += -openmp
LFLAGS += -openmp
else
DEFS += -DUSE_RANDOM
CFLAGS += -fopenmp -march=native -Wa,-q -Wvector-operation-performance
LFLAGS += -fopenmp -march=native -Wa,-q -Wvector-operation-performance
endif
endif

ifeq ($(INTEL_COMPILER), $(YES))
CFLAGS += -std=c++0x
LFLAGS += -std=c++0x
else
CFLAGS += -std=c++11
LFLAGS += -lstdc++.6
endif

export CFLAGS
export LFLAGS
export CC
export CL
export DEFS

all:
	cd gpc; make
	cd src; make

clean: 
	rm -f *~
	rm -f cv
	rm -f build/*.d
	rm -f build/*.o
	cd gpc; make clean
	cd src; make clean
