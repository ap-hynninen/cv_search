
OBJS = 	cv.o \
	cv_util.o \
	LM.o \
	Coord.o \
	GA.o \
	Gene.o \
	RW.o \
	Genome.o \
	cvGP.o

../cv : $(OBJS)
	$(CL) $(LFLAGS) -o ../cv $(OBJS) -L../gpc/lib -lgp
#	$(CL) $(LFLAGS) -o ../cv $(OBJS) -L../../gpc++0.5.2/lib -lgp

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

