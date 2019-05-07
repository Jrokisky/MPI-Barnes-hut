#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "octree.h"
#include "bucket.h"
#include "particle.h"

#define THETA 1.0
// Thanks newton.
#define GRAVITY 4.30091e-3

Octree * build_octree(OctantBucket *b, Space space)
{
    // Initialize tree.
    Octree * new_tree;
    new_tree = (Octree *) malloc(sizeof(Octree));
    // TODO: update size to use diag.
    new_tree->box_size = space.boundary_x - space.origin_x;
    new_tree->com_x = 0.0;
    new_tree->com_y = 0.0;
    new_tree->com_z = 0.0;
    new_tree->total_mass = 0.0;
    new_tree->num_leaves = 0;

    // Initialize tree.
    for (int i = 0; i < 8; i++) {
        (new_tree->children)[i] = NULL;
    }
    new_tree->value = NULL;

    // Initialize octant buckets.
    OctantBucket *buckets[8];
    for (int i = 0; i < 8; i++) {
        buckets[i] = NULL;
    }

    // Split particles from current bucket into its octant bucket.
    OctantBucket * traversal = b;
    while (traversal != NULL) {
        int octant = get_octant(*(traversal->value), space);
        add_to_bucket(&(buckets[octant]), traversal->value);
        traversal = traversal->next;
    }

    // Split current space into sub spaces.
    Space * sub_spaces[8];
    split_space(space, sub_spaces);

    for (int i = 0; i < 8; i++) {
        OctantBucket *curr_bucket = buckets[i];

        //Bucket has elements.
        if (curr_bucket != NULL) {

            // Bucket only has one element, so leaf.
            if (curr_bucket->next == NULL) {
                Octree * leaf = (Octree *) malloc(sizeof(Octree));
                leaf->value = curr_bucket->value;
                (new_tree->children)[i] = leaf;

                for (int i = 0; i < 8; i++) {
                    (leaf->children)[i] = NULL;
                }
                leaf->box_size = space.boundary_x - space.origin_x;
                (new_tree->num_leaves)++;
            }
            else {
                (new_tree->children)[i] = build_octree(curr_bucket, *(sub_spaces[i]));
                (new_tree->num_leaves) += ((new_tree->children)[i])->num_leaves;
            }
            free_bucket(curr_bucket);
        }
        free(sub_spaces[i]);
    }
    
    return new_tree;
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
        printf("\n");
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

void compute_center_of_mass(Octree *octree)
{
    if (octree == NULL) {
        return;
    }
    // Leaf
    else if (octree->value != NULL) {
        octree->total_mass = octree->value->mass;
        octree->com_x = octree->value->x;
        octree->com_y = octree->value->y;
        octree->com_z = octree->value->z;
        return;
    }
    
    // Compute total mass
    double tmp_mass = 0.0;
    double x_sum, y_sum, z_sum;
    x_sum = y_sum = z_sum = 0.0;
    for (int i = 0; i < 8; i++) {
        Octree * child = (octree->children)[i];
        if(child != NULL) {
            compute_center_of_mass(child);
            tmp_mass += clamp(child->total_mass); 
            x_sum += clamp(child->total_mass * child->com_x);
            y_sum += clamp(child->total_mass * child->com_y);
            z_sum += clamp(child->total_mass * child->com_z);
        }
    }
    octree->total_mass = tmp_mass;
    octree->com_x = clamp(x_sum / tmp_mass);
    octree->com_y = clamp(y_sum / tmp_mass);
    octree->com_z = clamp(z_sum / tmp_mass);
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

void compute_all_forces(Octree *traverser, Octree *whole_tree)
{
    if (traverser == NULL) {
        return;
    }
    else if (traverser->value != NULL) {
        compute_force(traverser->value, whole_tree);
    }
    else {
         for (int i = 0; i < 8; i++) {
            if ((traverser->children)[i] != NULL) {
                compute_all_forces((traverser->children)[i], whole_tree);
            }
        }
    }
}

void update_particle_position_and_velocity(Particle *p)
{
    double dt = 0.001;
    double acc_x = clamp(p->force_x / p->mass);
    double acc_y = clamp(p->force_y / p->mass);
    double acc_z = clamp(p->force_z / p->mass);

    // Compute middle of timestep velocity.
    double v_mid_x = p->vel_x + clamp( 0.5 * dt * acc_x);
    double v_mid_y = p->vel_y + clamp( 0.5 * dt * acc_y);
    double v_mid_z = p->vel_z + clamp( 0.5 * dt * acc_z);

    // Compute new position.
    double prev = p->x;
    p->x = p->x + clamp( dt * v_mid_x);
    p->y = p->y + clamp( dt * v_mid_y);
    p->z = p->z + clamp( dt * v_mid_z);

    // Compute second half of velocity.
    p->vel_x = v_mid_x + clamp( 0.5 * dt * acc_x );    
    p->vel_y = v_mid_y + clamp( 0.5 * dt * acc_y );    
    p->vel_z = v_mid_z + clamp( 0.5 * dt * acc_z );
}


void update_position_and_velocity(Octree *octree)
{
    if (octree == NULL) return;
    else if (octree->value != NULL) {
        Particle *p = octree->value;
        update_particle_position_and_velocity(p);
    }
    else {
        for (int i = 0; i < 8; i++) {
            Octree *child = (octree->children)[i];
            if (child != NULL) {
                update_position_and_velocity(child);
            }
        }
    }
}

void get_leaves(Octree *octree, Particle *ps, int * p_idx)
{
    if (octree == NULL) {
        return;
    }
    else if (octree->value != NULL) {
        memcpy(ps+(*p_idx), octree->value, sizeof(Particle));
        (*p_idx)++;
    }
    else {
         for (int i = 0; i < 8; i++) {
            if ((octree->children)[i] != NULL) {
                get_leaves((octree->children)[i], ps, p_idx);
            }
        }
    }
}
