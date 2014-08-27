
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

OBJS = cv.o cv_util.o LM.o Coord.o GA.o Pair.o

CFLAGS = -O3 
LFLAGS = -O3

ifeq ($(DEFS),USE_OPENMP)
ifeq ($(INTEL_COMPILER),1)
CFLAGS += -openmp
LFLAGS += -openmp
else
CFLAGS += -fopenmp -march=native -Wa,-q -Wvector-operation-performance
LFLAGS += -fopenmp -march=native -Wa,-q -Wvector-operation-performance
endif
endif

ifeq ($(INTEL_COMPILER),1)
CFLAGS += -std=c++0x
LFLAGS += -std=c++0x
else
CFLAGS += -std=c++11
LFLAGS += -lstdc++.6
endif

all: cv

# Link
cv : $(OBJS)
	$(CL) $(LFLAGS) -o cv $(OBJS)

clean: 
	rm -f *.o
	rm -f *~
	rm -f cv
	rm -f *.d

# Pull in dependencies that already exist
-include $(OBJS:.o=.d)

# Compile
%.o : %.cpp
	$(CC) -c $(CFLAGS) -D$(DEFS) $<
	$(CC) -MM $(CFLAGS) -D$(DEFS) $*.cpp > $*.d

