# Barnes-Hut Simulation using MPI

More information about the Barnes-Hut Simulation can be found [here](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation).

This implementation contains a serial and MPI version of the simulation and produces a gif of the simulated particle movement.

## Setup

The `gifmerge` executable needs to be compiled first. This can be done by running:

```
gcc gifmerge.c gifmerge.h -o gifmerge
```

## Running

The following scripts will generate a gif of the simulation in the `anim.gif` file

To run the serial version:

```
./run_serial_barnes_hut.sh $NUMBER_OF_PARTICLES $NUMBER_OF_SIMULATION_STEPS
```

To run the MPI version:
```
./run_mpi_barnes_hut.sh $NUMBER_OF_MPI_PROCESSES $NUMBER_OF_PARTICLES $NUMBER_OF_SIMULATION_STEPS
```
