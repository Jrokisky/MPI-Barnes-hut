#!/bin/bash
#Workflow script for particle sim
mydir=$PWD
rm anim.gif
rm timedat.0
g++ -g serial_barnes_hut.c particle.c octree.c -o serial-barnes-hut
./serial-barnes-hut $1 $2

#performn visualization
test=0;
while (($test <= $2-1))
do
echo "set xrange [0.0:150.0]
set title \"$1 particles at timestep $test\"
set yrange [0.0:150.0]
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

