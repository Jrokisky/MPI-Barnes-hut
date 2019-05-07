#include <stdlib.h>
#include <stdio.h>

#include "bucket.h"
#include "particle.h"

void add_to_bucket(OctantBucket **b, Particle *p)
{
    OctantBucket * new_item;
    new_item = (OctantBucket *) malloc(sizeof(OctantBucket));
    new_item->value = p;
    new_item->next = *b;
    *b = new_item;
}   

void print_bucket(OctantBucket *b)
{
    OctantBucket * curr = b;
    while (curr != NULL)
    {
        printf("[ ");
        print_particle(*(curr->value));
        printf(" ] => ");
        curr = curr->next;
    }
    printf("[ ]\n");
}

void free_bucket(OctantBucket *b)
{
    OctantBucket *curr = b;
    while (curr != NULL) {
        OctantBucket *tmp = curr->next;
        free(curr);
        curr = tmp;
    }
}

