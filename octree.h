#ifndef BARNES_HUT_OCTREE
#define BARNES_HUT_OCTREE

#include "particle.h"
#include "bucket.h"

typedef struct Octree {
    struct Octree *children[8];
    Particle *value;
    double total_mass;
    double com_x, com_y, com_z;
    double box_size;
    int num_leaves;
} Octree;

/**
 * Build an Octree from a set of particles.
 */
Octree * build_octree(OctantBucket *b, Space space);

void compute_center_of_mass(Octree *octree);

void compute_force(Particle *leaf, Octree *octree);

void compute_all_forces(Octree *traverser, Octree *octree);

void update_position_and_velocity(Octree *octree);

void update_particle_position_and_velocity(Particle *p);

void get_leaves(Octree *subtree, Particle *ps, int * p_idx);

void print_octree(Octree *octree, int indent);

void free_octree(Octree *octree);

#endif
