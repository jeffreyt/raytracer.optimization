#!/bin/bash

OUT="final_results.out"
optList="-O0 -O2 -O3 -O4"
T="1"

echo "testing started" > $OUT
echo "===================================" >> $OUT

for opt in $optList
	do
		echo "optimization $opt" >> $OUT
		gcc -DNUM_THREADS=$T -fopenmp -ffast-math $opt -c ezxml.c -o ezxml.o
		gcc -DNUM_THREADS=$T -fopenmp -ffast-math $opt -c graphics.cpp -o graphics.o
		g++ -DNUM_THREADS=$T -fopenmp -ffast-math $opt -c main.cpp -o main.o
		g++ -DNUM_THREADS=$T -fopenmp -ffast-math $opt ezxml.o graphics.o main.o -o melange
		for i in {1..5}
			do
				./melange >> $OUT
			done
		echo "===================================" >> $OUT
	done
