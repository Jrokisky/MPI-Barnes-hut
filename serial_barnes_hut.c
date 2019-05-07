#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "particle.h"
#include "bucket.h"
#include "octree.h"

#define WIDTH 100.0
#define LENGTH 100.0
#define HEIGHT 100.0

int main(int argc, char *argv[]) {
    srand(time(NULL));
    Space space = {WIDTH, LENGTH, HEIGHT, 0.0, 0.0, 0.0};
    OctantBucket * bucket = NULL;
    OctantBucket * head = NULL;
    int npart, t_step;

    if (argc < 2) {  
	    fprintf( stderr, "Usage: %s n timesteps\n", argv[0] ); 
    } 
    npart = atoi(argv[1]); 
    t_step = atoi(argv[2]);

    FILE *fp;
    fp = fopen("timedat.0", "a");

    Particle *particle_array = (Particle *) malloc(npart * sizeof(Particle));
    // Generate random particles.
    generate_random_particles(particle_array, space, npart);
    for (int i = 0; i < npart; i++) {
        print_particle(particle_array[i]);
        printf("\n");
        add_to_bucket(&bucket, &(particle_array[i]));
    }
    head = bucket;
    print_bucket(bucket);

    for (int i = t_step; i > 0; i--) {

        //    print_bucket(bucket);
        Octree * octree = build_octree(bucket, space);
        compute_center_of_mass(octree);
        compute_all_forces(octree, octree);
//        print_octree(octree, 0);
        update_position_and_velocity(octree);
        free_octree(octree);
    
        // Print particles in proper format.
        int p_id = 0;
        while (bucket != NULL) {
            Particle p = *(bucket->value);
            fprintf(fp,"%d %d %d %f %f %f\n",0, p_id, i, p.x, p.y, p.z);
            bucket = bucket->next;
            p_id++;
        }
        fprintf(fp, "\n\n");
        bucket = head;
    }
    
    fclose(fp);
    free(particle_array);

    // Free bucket and particles.
    free_bucket(bucket);
    return 0;
}
