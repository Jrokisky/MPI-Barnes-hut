#!/bin/bash

#Workflow script for particle sim
mydir=$PWD

if [[ -f "anim.gif" ]]; then
	rm anim.gif
fi

if [[ -f "timedat.0" ]]; then
	rm timedat.0
fi

#Compile executable
g++ -g serial_barnes_hut.c particle.c octree.c -o serial-barnes-hut
./serial-barnes-hut $1 $2

#Generate visualization plots for gif creation.
i=0;
while (($i <= $2-1))
do
echo "set xrange [0.0:100.0]
set title \"$1 particles at timestep $i\"
set yrange [0.0:100.0]
set zrange [0.0:100.0]
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

