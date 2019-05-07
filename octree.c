#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "octree.h"
#include "particle.h"

#define THETA 1.0
// Thanks newton.
#define GRAVITY 4.30091e-3

Octree * create_empty_octree(Space space)
{
    Octree *octree = (Octree *) malloc(sizeof(Octree));
    for (int i = 0; i < 8; i++) {
        (octree->children)[i] = NULL;
    }
    octree->value = NULL;
    octree->total_mass = 0.0;
    octree->com_x = 0.0;
    octree->com_y = 0.0;
    octree->com_z = 0.0;
    octree->box_size = space.boundary_x - space.origin_x;
    octree->num_leaves = 0;
    return octree;
}

void octree_insert(Octree *octree, Space space, Particle *p){
    bool done = false;
    // Originally used recursion, but trying to resolve performance issues.
    while (!done) 
    {

    update_center_of_mass(octree, p);
    // Check if we have any children.
    bool has_children = false;
    for (int i = 0; i < 8; i++) {
        if ((octree->children)[i] != NULL) {
            has_children = true;
            break;
        }
    }
    // Build subspaces
    Space sub_spaces[8];
    double mid_x = (space.origin_x + space.boundary_x)/2.0;
    
    double mid_y = (space.origin_y + space.boundary_y)/2.0;
    double mid_z = (space.origin_z + space.boundary_z)/2.0;

    // BOTTOM    TOP
    //  __ __   __ __
    // |2 |3 | |6 |7 |
    // |__|__| |__|__|
    // |0 |1 | |4 |5 |
    // |__|__| |__|__|
        
    for (int i = 0; i < 8; i++) {
        // Check if this octant is greater than the midpoint on different axis.
        bool gt_mid_x = i & (1 << 0);
        bool gt_mid_y = i & (1 << 1);
        bool gt_mid_z = i & (1 << 2);

        sub_spaces[i].boundary_x = gt_mid_x ? space.boundary_x : mid_x;
        sub_spaces[i].origin_x = gt_mid_x ? mid_x : space.origin_x;
        sub_spaces[i].boundary_y = gt_mid_y ? space.boundary_y : mid_y;
        sub_spaces[i].origin_y = gt_mid_y ? mid_y : space.origin_y;
        sub_spaces[i].boundary_z = gt_mid_z ? space.boundary_z : mid_z;
        sub_spaces[i].origin_z = gt_mid_z ? mid_z : space.origin_z;
    }

    if (octree->value != NULL) // LEAF.
    {
        // Stash the current particle and then clear it..
        Particle *curr_child = octree->value;
        octree->value = NULL;
        // Handle current child particle.        
        int curr_child_octant = get_octant(*curr_child, space);

        if ((octree->children)[curr_child_octant] == NULL) { 
            Octree * new_child  = create_empty_octree(sub_spaces[curr_child_octant]);
            (octree->children)[curr_child_octant] = new_child;
            update_center_of_mass(new_child, curr_child);
        }
        (octree->children)[curr_child_octant]->value = curr_child;

        // Insert new particle.
        int p_octant = get_octant(*p, space);

        // Old and new particle end up in different branches.
        if (p_octant != curr_child_octant) { 
            Octree * new_child  = create_empty_octree(sub_spaces[p_octant]);
            (octree->children)[p_octant] = new_child;
            update_center_of_mass(new_child, p);
            new_child->value = p;
            done = true;
        }
        else {
            // Same branch, so repeat.
            space = sub_spaces[p_octant];
            octree = (octree->children)[p_octant];
        }
    }
    else if(has_children) // MIDDLE LAYER.
    {
        // Insert new particle.
        int p_octant = get_octant(*p, space);
        if ((octree->children)[p_octant] == NULL) { 
            (octree->children)[p_octant] = create_empty_octree(sub_spaces[p_octant]);
        }
        space = sub_spaces[p_octant];
        octree = (octree->children)[p_octant];
    }
    else { // EMPTY LEAF (Root)
        octree->value = p;
        done = true;
    }

    }
}

void print_octree(Octree *octree, int indent) 
{
    for (int i = 0; i < indent; i++) {
        printf("    ");
    }

    printf("|-");
    if (octree == NULL) {
        printf("\n");
        return;
    }
    else if (octree->value == NULL) {
        printf("---");
        printf(" Mass: %f | CoM: (%f, %f, %f) | Box: %f | NumLeaves: %d\n",
                octree->total_mass,
                octree->com_x,
                octree->com_y,
                octree->com_z,
                octree->box_size,
                octree->num_leaves);
        for (int i = 0; i < 8; i++) {
            print_octree((octree->children)[i], indent+1);
        }
    }
    else {
        print_particle(*(octree->value));
    }
}    

void free_octree(Octree *octree)
{
    if (octree == NULL) {
        return;
    }
    for (int i = 0; i < 8; i++) {
        if ((octree->children)[i] != NULL) {
            free_octree((octree->children)[i]);
        }
    }
    free(octree);
}

void update_center_of_mass(Octree *octree, Particle *p)
{
    double total_mass = octree->total_mass + p->mass;
    octree->com_x = clamp( clamp(octree->total_mass * octree->com_x) + clamp(p->x * p->mass)) / total_mass;
    octree->com_y = clamp( clamp(octree->total_mass * octree->com_y) + clamp(p->y * p->mass)) / total_mass;
    octree->com_z = clamp( clamp(octree->total_mass * octree->com_z) + clamp(p->z * p->mass)) / total_mass;
}

void compute_force(Particle *leaf, Octree *octree)
{
    // When we get to a leaf.
    if (octree->value != NULL) {
        Particle *l = octree->value;
        double distance = clamp(compute_distance(leaf, l->x, l->y, l->z));

        leaf->force_x += clamp(GRAVITY * leaf->mass * l->mass * clamp((l->x - leaf->x) / clamp(pow(distance, 3.0))));
        leaf->force_y += clamp(GRAVITY * leaf->mass * l->mass * clamp((l->y - leaf->y) / clamp(pow(distance, 3.0))));
        leaf->force_z += clamp(GRAVITY * leaf->mass * l->mass * clamp((l->z - leaf->z) / clamp(pow(distance, 3.0))));
    }
    else {
        double distance = clamp(compute_distance(leaf, octree->com_x, octree->com_y, octree->com_z));

        // Use center of mass of octant.
        if ((octree->box_size / distance) < THETA) {
            leaf->force_x += clamp(GRAVITY * leaf->mass * octree->total_mass * clamp((octree->com_x - leaf->x) / clamp(pow(distance, 3.0))));
            leaf->force_y += clamp(GRAVITY * leaf->mass * octree->total_mass * clamp((octree->com_y - leaf->y) / clamp(pow(distance, 3.0))));
            leaf->force_z += clamp(GRAVITY * leaf->mass * octree->total_mass * clamp((octree->com_z - leaf->z) / clamp(pow(distance, 3.0))));
        }
        else {
            for (int i = 0; i < 8; i++) {
                if ((octree->children)[i] != NULL) {
                    compute_force(leaf, (octree->children)[i]);
                } 
            }
        }
    }
}

