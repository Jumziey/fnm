#ifndef __WAVEPACKET_PRIVATE__
#define __WAVEPACKET_PRIVATE__

#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>

/* Constants */

#define PISQR 9.86960440108935861883449099988 /* pi squared */

#define HALF_STEP 0
#define FULL_STEP 1

/* Macros */

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQUARE(z) ((z) * (z))

/* Structures */

typedef struct
{
  /* Flags indicating values to be calculated and printed out */
  bool yes, norm, energy, x_avg, y_avg, z_avg, sx, sy, sz, autoc, user;
  unsigned int wf;
} output_flags;

typedef struct
{
  /* Precomputed data needed for the application 
     of the kinetic energy operator */
  fftw_plan forward, backward;
  size_t n1, n1_local, n1_0, n2, n3;
  double complex *ek1, *ek2, *ek3, *ek1_2, *ek2_2, *ek3_2;
} kin_data;


typedef struct
{
  /* Precomputed data needed for the calculation of the energy */ 
  fftw_plan forward, backward;
  size_t n1, n1_local, n1_0, n2, n3;
  double *k1, *k2, *k3;
} ener_data;


/* Function declarations */

double energy (const parameters, const double, const ener_data, 
	       double * const, double complex *, double complex *);

void initialize_energy (const parameters, double complex *, 
			double complex *, ener_data *);

void initialize_kinetic (const parameters, double complex *, kin_data *);

void kinetic_operator (const kin_data, double complex *, const int);

void make_grid (parameters *);

void observe (const parameters, const output_flags, const double, 
	      const ener_data, double * const, const double complex * const, 
	      double complex *, double complex *, FILE *);

void potential_operator (const parameters, const double complex,
			 const double, double * const, double complex * const);

void print_header (const parameters, const output_flags, const char * const, 
		   FILE **);

void print_wf (const parameters, const double complex * const, 
	       const char * const);

void print_wf_bin (const parameters, const double complex * const, 
		   const char * const); 

void read_parameters (parameters * const, output_flags * const, 
		      unsigned int * const, unsigned int * const, 
		      const char * const, char * const, char * const, 
		      char * const);

#ifdef MPI
void distribute_parameters (parameters *, output_flags *, unsigned int *, 
		      unsigned int *);
#endif

#endif /* __WAVEPACKET_PRIVATE__ */

