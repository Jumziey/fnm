#
# Makefile for program "wavepacket"
#
# C.M. Dion, A. Hashemloo, and G. Rahali
# April 2013
#

CC = gcc

CFLAGS = -O3 -std=c99 

LDLIBS = -lfftw3 -lm 

freeEvoWave:	wavepacket.o freeEvo.o
	$(CC) -o freeEvoWave $^ $(LDLIBS)
	
potWellWave:	wavepacket.o potWell.o
	$(CC) -o potWellWave $^ $(LDLIBS)
	
squareBarrWave:	wavepacket.o squareBarr.o
	$(CC) -o squareBarrWave $^ $(LDLIBS)


