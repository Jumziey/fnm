#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h> //LAZY!
#include "wavepacket.h"

void
initialize_potential (const parameters params, const int argv, char ** const argc){}
void
initialize_user_observe (const parameters params, const int argc, char ** const argv){}

void
user_observe (const parameters params, const double t, const double complex * const psi){}

void
potential(const parameters p, const double t, double* const pot)
{
	int hwp,i;
	double v0 = 5.0e4, width=0.032;
		
	//half of well points
	hwp = (width/p.dx)/2;
	
	for(i=0; i<p.nx; i++)
			pot[i] = 0;
}

void
initialize_wf(const parameters p, const int argv, char** const argc, double complex *psi)
{
	int i;
	double k0 = 200;
	double x0 = 1.0/4.0;
	double s0 = 1.0/40.0;
	
	for(i=0; i<p.nx; i++)
		psi[i] = pow(1/(M_PI*s0*s0),0.25) * cexp(I*k0*p.x[i])*cexp(-pow((p.x[i]-x0),2)/(2*s0*s0));		
}
