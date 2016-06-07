CPP = g++
CC = gcc

NUM_THREADS = 2

CFLAGS= -DNUM_THREADS=$(NUM_THREADS) -O4

HEADERS = ezxml.h graphics.h melange.h

default: melange
	chmod +x ./compare_ppm.pl

melange: ezxml.o graphics.o main.o
	$(CPP) $(CFLAGS) ezxml.o graphics.o main.o -o melange

ezxml.o: ezxml.c $(HEADERS)
	$(CC) $(CFLAGS) -c ezxml.c -o ezxml.o
	
main.o: main.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c main.cpp -o main.o
	
graphics.o: graphics.cpp $(HEADER)
	$(CC) $(CFLAGS) -c graphics.cpp -o graphics.o

run: melange
	./melange

clean: 
	/bin/rm -f *.o melange rendered.ppm
