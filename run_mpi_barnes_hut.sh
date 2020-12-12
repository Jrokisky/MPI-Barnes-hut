#!/bin/bash

#Workflow script for particle sim
mydir=$PWD
rm anim.gif
rm timedat.0
mpicc -o mpi-barnes-hut mpi_barnes_hut.c octree.c particle.c -lm
mpirun -np $1 ./mpi-barnes-hut $2 $3

#performn visualization
i=0;
while (($i <= $3-1))
do
echo "set xrange [0.0:100.0]
set title \"$2 particles at timestep $i\"
set yrange [0.0:100.0]
set timestamp top offset 15 
set grid
set term gif size 1280,1280
set output '$i.gif'
plot \"timedat.0\" i $i u 4:5 pt 3  ps 1 t \"Node 0\";" >data_$i.gnu
gnuplot data_$i.gnu
let i=$i+1
done

#cleanup
rm -f *.gnu
ls *.gif | sort -nk1 | xargs ./gifmerge -10 -l0 >anim.vid
rm -f *.gif
mv anim.vid anim.gif

