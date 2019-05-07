/******************************************************************************
 * Implementation of the Barnes Hut Algorithm using MPI.
 * 
 * Author: Justin Rokisky
 * To compile: mpicc -o mpi-barnes-hut mpi_barnes_hut.c bucket.c octree.c particle.c -lm
 * To run: ./mpi_gen.sh NUM_PARTICLES NUM_TIMESTEPS
 * 
 * Design decisions:
 *  -> Currently only runs with 8 processes
 *  -> Each process builds a copy of the tree on its machine
 *  -> Each process independently takes one of the 8 subtrees and updates the
 *     position and velocity of the leafs in that tree.
 *  -> The tree is constructed from the bottom up. At each level, particles are
 *     split in the octants they will end up in, and then recursively passed 
 *     into the build_octree function.My original plan was to parallelize the
 *     tree construction but was unable to build a data structure that could 
 *     accomplish passing the tree pieces.
 *
 *****************************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "particle.h"
#include "octree.h"

#define WIDTH 100.0
#define LENGTH 100.0
#define HEIGHT 100.0

int main(int argc, char *argv[]) {
    int npart, t_step, rank, size;
    double      starttime,endtime;
    int debug = 1;

    // Init MPI.
    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Datatype mpi_particle_type;

    // Create MPI_Datatype for sending particles.
    MPI_Aint displacements_c[11] = {0, 4, 12, 20, 28, 36, 44, 52, 60, 68, 76};
    int block_lengths_c[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype dtypes_c[11] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct( 11, block_lengths_c, displacements_c, dtypes_c, &mpi_particle_type );
    MPI_Type_commit( &mpi_particle_type );

    // Default space.
    Space space = {WIDTH, LENGTH, HEIGHT, 0.0, 0.0, 0.0};

    starttime = MPI_Wtime(); 
    // Process user input
    if (argc < 2) {  
	    fprintf( stderr, "Usage: %s n timesteps\n", argv[0] ); 
        MPI_Abort( MPI_COMM_WORLD, 1 );
    } 
    npart = atoi(argv[1]); 
    t_step = atoi(argv[2]);

    Particle *particle_array = (Particle *) malloc(npart * sizeof(Particle));
    for (int i = 0; i < npart; i++) {
        // init empty.
        particle_array[i].id = 0;
        particle_array[i].x = 0.0;
        particle_array[i].y = 0.0;
        particle_array[i].z = 0.0;
        particle_array[i].vel_x = 0.0;
        particle_array[i].vel_y = 0.0;
        particle_array[i].vel_z = 0.0;
        particle_array[i].force_x = 0.0;
        particle_array[i].force_y = 0.0;
        particle_array[i].force_z = 0.0;
        particle_array[i].mass = 0.0;
    }
    // Process 0 generates and distributes random particles.
    if ( rank == 0 ) {
        srand(time(NULL));
        // Generate random particles.
        generate_random_particles(particle_array, space, npart);

        // Broadcast particle data.
        MPI_Bcast(particle_array, npart, mpi_particle_type, 0, MPI_COMM_WORLD);
    }
    else {
        // Get random particles from root node. 
        MPI_Bcast(particle_array, npart, mpi_particle_type, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = t_step; i > 0; i--) {
        // Each process builds the whole tree and computes centers of mass.
        Octree * octree = create_empty_octree(space);
        for (int j = 0; j < npart; j++) {
            octree_insert(octree, space, &(particle_array[j]));
        }

        int chunk = 0;
        int last_chunk = 0;
        // Split particles into chunks.
        if (npart / size == 0 || size == 1) {
            chunk = npart / size;
            last_chunk = chunk;
        }
        else {
            chunk = npart / (size - 1);
            last_chunk = npart % size;
        }

        // Each process processes its chunk.
        for (int idx = rank * chunk; idx < npart; idx++) {
            Particle *tmp_p = &(particle_array[idx]);
            compute_force(tmp_p, octree);
            update_particle_position_and_velocity(tmp_p);
        }
        // No longer need the octree.
        free_octree(octree);

        // Share this processes updated particle data with everyone else.
        for (int q = 0; q < size; q++) {
            // last chunk.
            if (size > 1 && q == (size-1)) chunk = last_chunk;
            Particle * updated_particles = (Particle *) malloc(chunk * sizeof(Particle));

            if (q == rank) {
                for (int f = 0; f < chunk; f++) {
                    int idx = rank * chunk + f;
                    memcpy(updated_particles + f, &(particle_array[idx]), sizeof(Particle));
                }       
            }

            // Perform the broadcast.
            MPI_Bcast( updated_particles, chunk, mpi_particle_type, q, MPI_COMM_WORLD );

            if (q != rank) {
                // Update particles.
                for (int k = 0; k < chunk; k++) {
                    // Get the index of the particle that needs updating.
                    int idx = updated_particles[k].id;
                    if (particle_array[idx].id == idx) {
                        particle_array[idx].x = updated_particles[k].x;
                        particle_array[idx].y = updated_particles[k].y;
                        particle_array[idx].z = updated_particles[k].z;
                        particle_array[idx].vel_x = updated_particles[k].vel_x;
                        particle_array[idx].vel_y = updated_particles[k].vel_y;
                        particle_array[idx].vel_z = updated_particles[k].vel_z;

                    }
                    else {
                        printf("WOT\n");
                    }
                }
            }
            free(updated_particles);
        }

        if (rank == 0) {
            FILE *fp;
            fp = fopen("timedat.0", "a");
            // Print particles in proper format.
            for (int step = 0; step < npart; step++) {
                Particle p = particle_array[step];
                fprintf(fp,"%d %d %d %f %f %f\n",0, p.id, i, p.x, p.y, p.z);
            }
            fprintf(fp, "\n\n");
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
   
    // Free Particles. 
    free(particle_array);

    endtime = MPI_Wtime() - starttime; 
    if(rank == 0) printf("Run Time: %f | Num Particles: %d | Num Steps: %d\n", endtime, npart, t_step);
    MPI_Finalize();
    return 0;
}
