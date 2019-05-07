#ifndef BARNES_HUT_OCTREE
#define BARNES_HUT_OCTREE

#include "particle.h"

typedef struct Octree {
    struct Octree *children[8];
    Particle *value;
    double total_mass;
    double com_x, com_y, com_z;
    double box_size;
    int num_leaves;
    Space space;
} Octree;

Octree * create_empty_octree(Space space);

void octree_insert(Octree *octree, Space space, Particle *p);

void update_center_of_mass(Octree *octree, Particle *p);

void compute_force(Particle *leaf, Octree *octree);

void print_octree(Octree *octree, int indent);

void free_octree(Octree *octree);

#endif
