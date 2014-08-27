
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

SRC = cv.cpp cv_util.cpp LM.cpp Coord.cpp GA.cpp Pair.cpp

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
LFLAGS += -std=c++0x
else
CFLAGS += -std=c++11
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

cv.o: /usr/include/sys/time.h /usr/include/sys/cdefs.h
cv.o: /usr/include/sys/_symbol_aliasing.h
cv.o: /usr/include/sys/_posix_availability.h /usr/include/sys/_types.h
cv.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
cv.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
cv.o: /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h
cv.o: /usr/include/sys/_types/_timeval.h /usr/include/sys/_types/_time_t.h
cv.o: /usr/include/sys/_types/_suseconds_t.h
cv.o: /usr/include/sys/_types/_fd_setsize.h /usr/include/sys/_types/_fd_set.h
cv.o: /usr/include/sys/_types/_fd_clr.h /usr/include/sys/_types/_fd_isset.h
cv.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
cv.o: /usr/include/time.h /usr/include/_types.h /usr/include/_structs.h
cv.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_null.h
cv.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_size_t.h
cv.o: /usr/include/sys/_select.h /usr/include/stdio.h
cv.o: /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h
cv.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
cv.o: /usr/include/secure/_common.h cv_util.h LM.hpp Coord.hpp struct.h
cv.o: GA.hpp Pair.hpp
cv_util.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
cv_util.o: /usr/include/sys/_symbol_aliasing.h
cv_util.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
cv_util.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
cv_util.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
cv_util.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
cv_util.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
cv_util.o: /usr/include/sys/_types/_off_t.h
cv_util.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
cv_util.o: /usr/include/secure/_common.h /usr/include/stdlib.h
cv_util.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
cv_util.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
cv_util.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
cv_util.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
cv_util.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
cv_util.o: /usr/include/sys/_types/_sigaltstack.h
cv_util.o: /usr/include/sys/_types/_ucontext.h
cv_util.o: /usr/include/sys/_types/_pthread_attr_t.h
cv_util.o: /usr/include/sys/_types/_sigset_t.h
cv_util.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
cv_util.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
cv_util.o: /usr/include/sys/_types/_int16_t.h
cv_util.o: /usr/include/sys/_types/_int32_t.h
cv_util.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
cv_util.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
cv_util.o: /usr/include/_types/_uint64_t.h
cv_util.o: /usr/include/sys/_types/_intptr_t.h
cv_util.o: /usr/include/sys/_types/_uintptr_t.h
cv_util.o: /usr/include/_types/_intmax_t.h /usr/include/_types/_uintmax_t.h
cv_util.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
cv_util.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
cv_util.o: /usr/include/libkern/_OSByteOrder.h
cv_util.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
cv_util.o: /usr/include/sys/_types/_ct_rune_t.h
cv_util.o: /usr/include/sys/_types/_rune_t.h
cv_util.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
cv_util.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
cv_util.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
cv_util.o: /usr/include/math.h
LM.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
LM.o: /usr/include/sys/_symbol_aliasing.h
LM.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
LM.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
LM.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
LM.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
LM.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
LM.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
LM.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
LM.o: /usr/include/stdlib.h /usr/include/sys/wait.h
LM.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
LM.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
LM.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
LM.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
LM.o: /usr/include/mach/i386/_structs.h
LM.o: /usr/include/sys/_types/_sigaltstack.h
LM.o: /usr/include/sys/_types/_ucontext.h
LM.o: /usr/include/sys/_types/_pthread_attr_t.h
LM.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
LM.o: /usr/include/sys/resource.h /usr/include/stdint.h
LM.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
LM.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
LM.o: /usr/include/_types/_uint8_t.h /usr/include/_types/_uint16_t.h
LM.o: /usr/include/_types/_uint32_t.h /usr/include/_types/_uint64_t.h
LM.o: /usr/include/sys/_types/_intptr_t.h
LM.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
LM.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
LM.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
LM.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
LM.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
LM.o: /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h
LM.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
LM.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
LM.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
LM.o: /usr/include/math.h LM.hpp
Coord.o: /usr/include/stdlib.h /usr/include/Availability.h
Coord.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
Coord.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
Coord.o: /usr/include/sys/_symbol_aliasing.h
Coord.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Coord.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
Coord.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Coord.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Coord.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Coord.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Coord.o: /usr/include/mach/i386/_structs.h
Coord.o: /usr/include/sys/_types/_sigaltstack.h
Coord.o: /usr/include/sys/_types/_ucontext.h
Coord.o: /usr/include/sys/_types/_pthread_attr_t.h
Coord.o: /usr/include/sys/_types/_sigset_t.h
Coord.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_uid_t.h
Coord.o: /usr/include/sys/resource.h /usr/include/stdint.h
Coord.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
Coord.o: /usr/include/sys/_types/_int32_t.h
Coord.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
Coord.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
Coord.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
Coord.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
Coord.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
Coord.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Coord.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Coord.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Coord.o: /usr/include/sys/_types/_ct_rune_t.h
Coord.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Coord.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
Coord.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
Coord.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Coord.o: Coord.hpp struct.h
GA.o: cv_util.h LM.hpp GA.hpp Coord.hpp struct.h Pair.hpp
Pair.o: Pair.hpp struct.h
