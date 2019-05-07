#!/bin/bash
#Workflow script for particle sim
mydir=$PWD
rm anim.gif
rm -f timedat.*

#execute code
mpirun -np 4 ./nbodypipe $1 $2

#performn visualization
test=0;
while (($test <= $2-1))
do
echo "set xrange [0:1]
set title \"$1 particles at timestep $test\"
set yrange [0:1]
set timestamp top offset 15 
set grid
set term gif size 320,480
set output '$test.gif'
plot \"timedat.0\" i $test u 4:5 pt 3  ps 1 t \"Node 0\", \"timedat.1\" i $test u 4:5 pt 3  ps 1 t \"Node 1\",\"timedat.2\" i $test u 4:5 pt 3 ps 1 t \"Node 2\", \"timedat.3\" i $test u 4:5 pt 3  ps 1 t \"Node 3\";" >data_$test.gnu
gnuplot data_$test.gnu
let test=$test+1
done


#cleanup
rm -f *.gnu
ls *.gif | sort -nk1 | xargs ./gifmerge -10 -l0 >anim.vid
rm -f *.gif
mv anim.vid anim.gif

