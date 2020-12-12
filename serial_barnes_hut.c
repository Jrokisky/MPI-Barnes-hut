// Serial implementation of the MPI Barnes Hut algorithm.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <float.h>

#include "particle.h"
#include "octree.h"

#define WIDTH 99999.0
#define LENGTH 99999.0
#define HEIGHT 99999.0

#define ORIGIN_X -99999.0
#define ORIGIN_Y -99999.0
#define ORIGIN_Z -99999.0

int main(int argc, char *argv[]) {
    int npart, t_step, rank, size;

    int debug = 1;

    // Default space.
    Space space = {WIDTH, LENGTH, HEIGHT, ORIGIN_X, ORIGIN_Y, ORIGIN_Z}; 

    // Process user input
    if (argc < 2) {  
	    fprintf( stderr, "Usage: %s n timesteps\n", argv[0] ); 
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
    srand(time(NULL));
    // Generate random particles.
    generate_random_particles(particle_array, space, npart);

    for (int i = t_step; i > 0; i--) {
        // Each process builds the whole tree and computes centers of mass.
        Octree * octree = create_empty_octree(space);
        for (int j = 0; j < npart; j++) {
            if (in_space(space,particle_array[j])) {
                octree_insert(octree, space, &(particle_array[j]));
            }
        }

        // Each process processes its chunk.
        for (int idx = 0; idx < npart; idx++) {
            Particle *tmp_p = &(particle_array[idx]);
            if (in_space(space, particle_array[idx])) {
                compute_force(tmp_p, octree);
                update_particle_position_and_velocity(tmp_p);
            }
        }
        // No longer need the octree.
        free_octree(octree);

        FILE *fp;
        fp = fopen("timedat.0", "a");
        // Print particles in proper format.
        for (int step = 0; step < npart; step++) {
            Particle p = particle_array[step];
            if (in_space(space, p)) {
                fprintf(fp,"%d %d %d %f %f %f\n",0, p.id, i, p.x, p.y, p.z);
            }
        }
        fprintf(fp, "\n\n");
        fclose(fp);
        printf("STEP: %d completed\n", i);
    }
   
    // Free Particles. 
    free(particle_array);
    return 0;
}
