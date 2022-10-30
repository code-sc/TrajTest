ODIR	 = obj
_OBJ	 = g_trajtest.o
OBJS	 = $(patsubst %,$(ODIR)/%,$(_OBJ))

IDIR	 = ../include
_DEP	 =
DEPS	 = $(patsubst %,$(IDIR)/%,$(_DEP))

LIBS	 = -lgromacs_d -lm -lfftw3 -pthread
EXE	 = g_trajtest

CC	 = mpicxx
CFLAGS	 = -O3 -finline-functions -g -Wall -pthread -std=c++1z -I/home/scj/Softwares/gromacs-double/include/ -DGMX_DOUBLE=1 -fopenmp
LDFLAGS	 = -g -Wall -std=c++1z -L/home/scj/Softwares/gromacs-double/lib/

$(ODIR)/%.o	: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXE)	: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean	:
	rm $(OBJS) $(EXE)
