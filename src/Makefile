
OBJS = 	cv.o \
	cv_util.o \
	LM.o \
	Coord.o \
	GA.o \
	Pair.o \
	RW.o \
	Genome.o

../cv : $(OBJS)
	$(CL) $(LFLAGS) -o ../cv $(OBJS)

clean: 
	rm -f *~
	rm -f *.o
	rm -f *.d

# Pull in dependencies that already exist
-include $(OBJS:.o=.d)

# Compile
%.o : %.cpp
	$(CC) -c $(CFLAGS) $(DEFS) $<
	$(CC) -MM $(CFLAGS) $(DEFS) $*.cpp > $*.d

