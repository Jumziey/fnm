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
	
potWellWavek0: wavepacket.o potWellk0.o
	$(CC) -o potWellWavek0 $^ $(LDLIBS)
	
squareBarrWave:	wavepacket.o squareBarr.o
	$(CC) -o squareBarrWave $^ $(LDLIBS)
	
clean:
	rm *.o squareBarrWave potWellWavek0 potWellWave freeEvoWave


