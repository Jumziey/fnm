#ifndef __WAVEPACKET__
#define __WAVEPACKET__

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

/* Structures */

typedef struct
{
  /* Parameters and grid */
  int size, rank;
  size_t nx, ny, nz, n, nx_local, nx0, n_local;
  double x_min, y_min, z_min, x_max, y_max, z_max, dx, dy, dz;
  double *x, *y, *z, *x2, *y2, *z2;
  double mass, dt, hbar;
  bool is1D;
} parameters;


/* Function declarations */

void abort ();

void distribute_wf (const parameters, double complex * const,
		    double complex * const);

double expectation1D (const parameters, const int, const double * const, 
		      const double complex * const);

void initialize_potential (const parameters, const int, char ** const);

void initialize_user_observe (const parameters, const int, char ** const);

void initialize_wf (const parameters, const int, char ** const, 
		    double complex *);

double complex integrate3D (const parameters, const double complex * const,
			    const double complex * const);

double norm (const parameters, const double complex * const);

void potential (const parameters, const double, double * const);

void renormalize (const parameters, double complex *);

void user_observe (const parameters, const double, 
		   const double complex * const);

#endif /* __WAVEPACKET__ */

