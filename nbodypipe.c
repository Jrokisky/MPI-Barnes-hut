/* Nbody simulation using MPI
*  Tyler Simon
*/ 


#include "mpi.h" 
#include <unistd.h>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h> 
#include <time.h>
#include <string.h>

#include "structs.h"

#define MAX_PARTICLES 8000 
#define MAX_P          128 
#define DEBUG 1
 
/* 
extern double drand48(); 
*/ 
 
/* Pipeline version of the algorithm... */ 
/* we really need the velocities as well... */ 
 
 
void InitParticles( int, Particle[], ParticleV [], int ); 
double ComputeForces( Particle [], Particle [], ParticleV [], int ); 
double ComputeNewPos( Particle [], ParticleV [], int, double, MPI_Comm ); 
 
int main( int argc, char *argv[] ) 
{ 
    Particle  particles[MAX_PARTICLES];   /* Particles on ALL nodes */ 
    ParticleV pv[MAX_PARTICLES];          /* Particle velocity */ 
    Particle  sendbuf[MAX_PARTICLES],     /* Pipeline buffers */ 
	recvbuf[MAX_PARTICLES]; 
    MPI_Request request[2]; 
    int         counts[MAX_P],              /* Number on each processor */ 
 	        displs[MAX_P];              /* Offsets into particles */ 
    int         rank, size, npart, i, j, 
	        offset;                     /* location of local particles */ 
    int         totpart;                    /* total number of particles */ 
    MPI_Datatype particletype; 
    double      sim_t;                      /* Simulation time */ 
    double      starttime,endtime;                       /* Computation time */ 
    int         pipe, left, right, periodic; 
    int 	t_step; 
    MPI_Comm    commring; 
    MPI_Status  statuses[2]; 
 
    MPI_Init( &argc, &argv ); 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size ); 
 
/* Get the best ring in the topology */ 
    periodic = 1; 
    MPI_Cart_create( MPI_COMM_WORLD, 1, &size, &periodic, 1, &commring ); 
    MPI_Cart_shift( commring, 0, 1, &left, &right ); 
 
/* Everyone COULD have a different size ... */ 
    if (argc < 2) {  
	fprintf( stderr, "Usage: %s n timesteps\n", argv[0] ); 
	MPI_Abort( MPI_COMM_WORLD, 1 ); 
    } 
    npart = atoi(argv[1]) / size; 
    t_step = atoi(argv[2]);
     
 
    if (npart * size > MAX_PARTICLES) { 
	fprintf( stderr, "%d is too many; max is %d\n",  
		 npart*size, MAX_PARTICLES ); 
	MPI_Abort( MPI_COMM_WORLD, 1 ); 
    } 
 
    MPI_Type_contiguous( 4, MPI_DOUBLE, &particletype ); 
    MPI_Type_commit( &particletype ); 
 
/* Get the sizes and displacements */ 
    MPI_Allgather( &npart, 1, MPI_INT, counts, 1, MPI_INT, commring ); 
    displs[0] = 0; 
    for (i=1; i<size; i++)  
	displs[i] = displs[i-1] + counts[i-1]; 
    totpart = displs[size-1] + counts[size-1]; 
 
/* Generate the initial values */ 
    srand48(time(NULL)*rank);
    InitParticles( rank,particles, pv, npart); 
    offset = displs[rank]; 
    FILE *fp;
    char wfile[255]; 
    starttime = MPI_Wtime(); 
    sim_t = 0.0; 

if(DEBUG){
sprintf(wfile,"timedat.%d", rank);
fp=fopen(wfile, "w+");
	}

    while (t_step--) { 
	double max_f, max_f_seg; 
     
	/* Load the initial sendbuffer */ 
	memcpy( sendbuf, particles, npart * sizeof(Particle) ); 
	max_f = 0.0; 
	for (pipe=0; pipe<size; pipe++) { 
	    if (pipe != size-1) { 
		MPI_Isend( sendbuf, npart, particletype, right, pipe,  
			   commring, &request[0] ); 
		MPI_Irecv( recvbuf, npart, particletype, left,  pipe,  
			   commring, &request[1] ); 
	    } 
	    /* Compute forces (2D only) */ 
	    max_f_seg = ComputeForces( particles, sendbuf, pv, npart ); 
	    if (max_f_seg > max_f) max_f = max_f_seg; 
	    /* Push pipe */ 
	    if (pipe != size-1)  
		MPI_Waitall( 2, request, statuses ); 
	    memcpy( sendbuf, recvbuf, counts[pipe] * sizeof(Particle) ); 
	} 
	/* Once we have the forces, we compute the changes in position */ 
	sim_t += ComputeNewPos( particles, pv, npart, max_f, commring ); 
if(DEBUG)
{
	for (i=1; i<npart; i++)
	  {
fprintf(fp,"%d %d %d %f %f %f\n",rank, i, t_step, particles[i].x, particles[i].y, particles[i].z);
	  printf("%d %d %d %f %f %f\n",rank, i, t_step, particles[i].x, particles[i].y, particles[i].z);
	  }
	  fprintf(fp,"\n\n");
	}//end debug

	 
 
	/* We could do graphics here (move particles on the display) */ 
	
    } 
if(DEBUG){
	fclose(fp);
}
    endtime = MPI_Wtime() - starttime; 
/*
    if (rank == 0) { 
	for (i=1; i<npart; i++)
	  {
	  fprintf(fp,"%d %d %d %f %f %f\n",rank, i, t_step, particles[i].x, particles[i].y, particles[i].z);
	  printf("%f %f %f\n",particles[i].x, particles[i].y, particles[i].z);
	  }
	
        printf("# rank particle no. timestep xpos ypos zpos\n");
	printf( "#Computed %d particles in %f seconds\n", totpart, endtime ); 
    } 
*/
    MPI_Finalize(); 
    return 0; 
} //end main
 
void InitParticles(int myrank, Particle particles[], ParticleV pv[], int npart ) 
{ 
    int i; 
    for (i=0; i<npart; i++) { 
	particles[i].x	  = drand48(); 
	particles[i].y	  = drand48(); 
	particles[i].z	  = drand48(); 

/* all of node 0 particles are at one spot and have more mass*/
	if(myrank==0)
	{
	particles[5].mass=10.0;
 	particles[5].x = 0.75;		
 	particles[5].y = 0.75;		
	particles[5].z	  = drand48(); 
	}
	else
	particles[i].mass=drand48();
	pv[i].xold	  = particles[i].x; 
	pv[i].yold	  = particles[i].y; 
	pv[i].zold	  = particles[i].z; 
	pv[i].fx	  = 0; 
	pv[i].fy	  = 0; 
	pv[i].fz	  = 0; 
    } 
} 
 
double ComputeForces( Particle myparticles[], Particle others[],  
		      ParticleV pv[], int npart ) 
{ 
  double max_f, rmin; 
  int i, j; 
 
  max_f = 0.0; 
  for (i=0; i<npart; i++) { 
    double xi, yi, zi, mi, rx, ry, rz, mj, r, fx, fy, fz; 
    rmin = 1; 
    xi   = myparticles[i].x; 
    yi   = myparticles[i].y; 
    zi   = myparticles[i].z; //velocity 10/17
    fx   = 0.0; 
    fy   = 0.0; 
    fz   = 0.0; //velocity 10/17 
    for (j=0; j<npart; j++) { 
      rx = xi - others[j].x; 
      ry = yi - others[j].y; 
      rz = zi - others[j].z; //velocity
      mj = others[j].mass; 
      r  = sqrt(rx * rx  + ry * ry); 
      /* ignore overlap and same particle */ 
      if (r == 0.0) continue; 
      if (r < rmin) rmin = r; 
      /* compute forces */ 
      r  = r * sqrt(r); 
      fx -= mj * rx / r; 
      fy -= mj * ry / r; 
      fz -= mj * rz / r;  //velocity
    } 
    pv[i].fx += fx; 
    pv[i].fy += fy; 
    pv[i].fz += fz; //velocity
    /* Compute a rough estimate of (1/m)|df / dx| */ 
    fx		      = sqrt(fx*fx + fy*fy)/rmin; 
    //fx		      = sqrt(fx*fx + fy*fy+ fz*fz)/rmin; //velocity
    if (fx > max_f) max_f = fx; 
	//printf("(%d) fx=%f\n",  i,fx);
  } 
  return max_f; 
} 
 
double ComputeNewPos( Particle particles[], ParticleV pv[], int npart,  
		      double max_f, MPI_Comm commring ) 
{ 
  int i; 
  double      a0, a1, a2; 
  static      double dt_old = 0.001, dt = 0.001; 
  double      dt_est, new_dt, dt_new; 
 
  /* integration is a0 * x^+ + a1 * x + a2 * x^- = f / m */ 
  a0	 = 2.0 / (dt * (dt + dt_old)); 
  a2	 = 2.0 / (dt_old * (dt + dt_old)); 
  a1	 = -(a0 + a2);      /* also -2/(dt*dt_old) */ 
 
  for (i=0; i<npart; i++) { 
    double xi, yi, zi; //velocity

    /* Very, very simple leapfrog time integration.  We use a variable  
       step version to simplify time-step control. */ 
    xi	           = particles[i].x; 
    yi	           = particles[i].y; 
    zi	           = particles[i].z; //velocity
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0; 
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0; 
    particles[i].z = (pv[i].fz - a1 * zi - a2 * pv[i].zold) / a0; //velocity
    pv[i].xold     = xi; 
    pv[i].yold     = yi; 
    pv[i].zold     = zi; 
    pv[i].fx       = 0; 
    pv[i].fy       = 0; 
    pv[i].fz       = 0; //velocity
  } 
 
  /* Recompute a time step. Stability criteria is roughly  
     2/sqrt(1/m |df/dx|) >= dt.  We leave a little room */ 
  dt_est = 1.0/sqrt(max_f); 
  /* Set a minimum: */ 
  if (dt_est < 1.0e-6) dt_est = 1.0e-6; 
  MPI_Allreduce( &dt_est, &dt_new, 1, MPI_DOUBLE, MPI_MIN, commring ); 
  /* Modify time step */ 
  if (dt_new < dt) { 
    dt_old = dt; 
    dt     = dt_new; 
  } 
  else if (dt_new > 4.0 * dt) { 
    dt_old = dt; 
    dt    *= 2.0; 
  } 
 
  return dt_old; 
} 

