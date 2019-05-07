#!/bin/bash
#Workflow script for particle sim
mydir=$PWD
rm anim.gif
rm timedat.0
mpicc -o mpi-barnes-hut mpi_barnes_hut.c octree.c particle.c -lm
mpirun -np $1 ./mpi-barnes-hut $2 $3

#performn visualization
test=0;
while (($test <= $3-1))
do
echo "set xrange [-50.0:100.0]
set title \"$1 particles at timestep $test\"
set yrange [-50.0:100.0]
set timestamp top offset 15 
set grid
set term gif size 1280,1280
set output '$test.gif'
plot \"timedat.0\" i $test u 4:5 pt 3  ps 1 t \"Node 0\";" >data_$test.gnu
gnuplot data_$test.gnu
let test=$test+1
done

#cleanup
rm -f *.gnu
ls *.gif | sort -nk1 | xargs ./gifmerge -10 -l0 >anim.vid
rm -f *.gif
mv anim.vid anim.gif

