# Comments in a Makefile starts with #

# Part 1. Definitins of tools flags etc
CFLAGS=-O2
LOADLIBES = -lgsl -lgslcblas -lm

# Part 2. Definitions of macros representing files
OBJECTS  = phonons.o frequencies.o

ALL_SOURCE = phonons.c frequencies.c

# Part 3. Build rules

phonons: $(OBJECTS)

# Part 5. Administrative targets

clean:
	$(RM) *~ *.o

print: $(ALL_SOURCE) Makefile
	lpr $?
	touch print

