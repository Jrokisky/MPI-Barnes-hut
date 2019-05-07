#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>

#include "particle.h"

double drand_custom (double low, double high)
{
    double range = high - low;
    return ((double)rand() * range) / (double)RAND_MAX + low;
}

double clamp(double in) 
{
    if ( in >= 0.0) {
        double tmp = in < DBL_MIN ? DBL_MIN : in;
        return tmp > DBL_MAX ? DBL_MAX : tmp;
    }
    else {
        double tmp = in < -DBL_MAX ? -DBL_MAX : in;
        return tmp > -DBL_MIN ? -DBL_MIN : tmp;
    }
}
int get_octant(Particle particle, Space space)
{
    int octant = 0;
    if (particle.x > ((space.origin_x + space.boundary_x)/2.0)) {
        octant |= (1 << 0);
    }

    if (particle.y > ((space.origin_y + space.boundary_y)/2.0)) {
        octant |= (1 << 1);
    }

    if (particle.z > ((space.origin_z + space.boundary_z)/2.0)) {
        octant |= (1 << 2);
    }
    return octant;
}

void split_space(Space space, Space *sub_spaces[])
{
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

        Space *octant = (Space *) malloc(sizeof(Space));
        octant->boundary_x = gt_mid_x ? space.boundary_x : mid_x;
        octant->origin_x = gt_mid_x ? mid_x : space.origin_x;
        octant->boundary_y = gt_mid_y ? space.boundary_y : mid_y;
        octant->origin_y = gt_mid_y ? mid_y : space.origin_y;
        octant->boundary_z = gt_mid_z ? space.boundary_z : mid_z;
        octant->origin_z = gt_mid_z ? mid_z : space.origin_z;
        sub_spaces[i] = octant;
    }
}

void generate_random_particles(Particle * p, Space space, int count) 
{
    for (int i = 0; i < count; i++) {
        double x = drand_custom(0.0, space.boundary_x);
        double y = drand_custom(0.0, space.boundary_y);
        double z = drand_custom(0.0, space.boundary_z);

        p[i].x = x;
        p[i].y = y;
        p[i].z = z;
        p[i].vel_x = 0.0;//drand_custom(-100.0, 100.0);
        p[i].vel_y = 0.0;//drand_custom(-100.0, 100.0);
        p[i].vel_z = 0.0;//drand_custom(-100.0, 100.0);
        p[i].mass = 10000.0;//drand_custom(10000.0, 100000.0);
        p[i].id = i;
    }
}

void print_particle(Particle p)
{
    printf("id: %d | location: (%f, %f, %f) | mass: %f | force: (%f, %f, %f)", p.id, p.x, p.y, p.z, p.mass, p.force_x, p.force_y, p.force_z);
}

void print_space(Space s)
{
    printf("origin: (%f, %f, %f) | boundary: (%f, %f, %f)", s.origin_x, s.origin_y, s.origin_z, s.boundary_x, s.boundary_y, s.boundary_z);
}

double compute_distance(Particle *particle, double com_x, double com_y, double com_z)
{
    return sqrt( clamp(pow(particle->x - com_x, 2)) + clamp(pow(particle->y - com_y, 2)) + clamp(pow(particle->z - com_z, 2)));
}



