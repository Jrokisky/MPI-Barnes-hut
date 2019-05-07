#ifndef BARNES_HUT_PARTICLES
#define BARNES_HUT_PARTICLES

typedef struct Particle { 
    int id;
    double x, y, z; 
    double vel_x, vel_y, vel_z;
    double force_x, force_y, force_z;
    double mass; 
} Particle; 

typedef struct Space {
    double boundary_x, boundary_y, boundary_z;
    double origin_x, origin_y, origin_z;
} Space;

double drand_custom(double low, double high);
double clamp(double in);
/**
 * Get the octant that the given particle belongs in.
 */
int get_octant(Particle particle, Space space);

/**
 * Split the given space into its octants.
 */
void split_space(Space space, Space *sub_spaces[]);

/**
 * Generate a random particle in the given space.
 */
void generate_random_particles(Particle * p, Space s, int count);

/**
 * Print a particle's information.
 */
void print_particle(Particle particle);

void print_space(Space s);

double compute_distance(Particle *particle, double com_x, double com_y, double com_z);

#endif
