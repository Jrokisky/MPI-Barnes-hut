#!/bin/bash

mpicc -o mpi-barnes-hut mpi_barnes_hut.c octree.c particle.c -lm
rm results.txt
for NUM_PROC in 1 2 4 8
do
    for NUM_PART in 100 200 400 800
	do
        for NUM_STEPS in 1000  
	    do
		    mpirun -np $NUM_PROC ./mpi-barnes-hut $NUM_PART $NUM_STEPS >> results.txt
		echo >> results.txt  # newline
        done
	done
done
echo >> results.txt
