#
# Makefile for program "wavepacket"
#
# C.M. Dion, A. Hashemloo, and G. Rahali
# April 2013
#

CC = gcc

CFLAGS = -O3 -std=c99 

LDLIBS = -lfftw3 -lm 

wavepacket:	wavepacket.o user_defined.o
	$(CC) -o wavepacket $^ $(LDLIBS)

.PHONY: clean
clean:	
	rm -f *.o iontrap3D co
