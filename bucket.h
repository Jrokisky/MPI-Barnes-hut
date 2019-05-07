#ifndef BARNES_HUT_BUCKET
#define BARNES_HUT_BUCKET

#include <stdio.h>

#include "particle.h"

// Linked list of Particles.
typedef struct OctantBucket {
    Particle *value;
    struct OctantBucket *next;
} OctantBucket;

/**
 * Add a Particle to appropriate bucket.
 */
void add_to_bucket(OctantBucket **b, Particle *p);

void print_bucket(OctantBucket *b);

void free_bucket(OctantBucket *b);

#endif
