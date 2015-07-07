#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_math.h> //LAZY!
#include "wavepacket.h"

double k0;

void
initialize_potential (const parameters params, const int argv, char ** const argc){}

void
initialize_user_observe (const parameters params, const int argc, char ** const argv)
{
	FILE* fp = fopen("potWellk0.dat","w+");
	fclose(fp);
}

void
user_observe(const parameters p, const double t, const double complex * const psi)
{
	int hwp,i;
	double width=0.032, T=0, R=0, dx;
	FILE* fp;
		
	//half of well points
	hwp = ceil((width/p.dx)/2);
	dx = (p.x_max-p.x_min)/p.nx;
	
	R += pow(cabs(psi[0]),2)*0.5;
	for(i=1; i<(p.nx-1); i++)
		if(i<ceil(p.nx/2)-hwp)
			R += pow(cabs(psi[i]),2);
		else if(i>ceil(p.nx/2)+hwp)
			T += pow(cabs(psi[i]),2);
	T += pow(cabs(psi[p.nx]),2)*0.5;
	
	T *= dx;
	R *= dx;
	
	if((fabs(t - 0.0011))<1e-12) {
		//Write down data
		fp = fopen("potWellk0.dat", "a");
		fprintf(fp,"%g %g %g\n", k0, T, R);
		fclose(fp);
		fprintf(stderr,"####################################WRITING\n");
	}
}

void
potential(const parameters p, const double t, double* const pot)
{
	int hwp,i;
	double v0 = -5.0e4, width=0.032;
		
	//half of well points
	hwp = ceil((width/p.dx)/2);
	
	for(i=0; i<p.nx; i++)
		if(i>ceil(p.nx/2)-hwp && i<ceil(p.nx/2)+hwp)
			pot[i] = v0;
		else
			pot[i] = 0;
}

void
initialize_wf(const parameters p, const int argc, char** const argv, double complex *psi)
{
	int i, opt;
	char* check;
	FILE* fp;

	double x0 = 1.0/4.0;
	double s0 = 1.0/40.0;
	
	
		while((opt=getopt(argc,argv,"k:")) != -1) {
		switch(opt) {
			case 'k':
				optarg++;//Removing equal sign
				k0 = strtol(optarg,&check,10);
				if(*check != '\0')
					fprintf(stderr,"only take integers as argument to k\n");
				break;
			case '?':
				continue;
			default:
				fprintf(stderr, "we shouldnt be here\n");
		}
	}
	
	for(i=0; i<p.nx; i++)
		psi[i] = pow(1/(M_PI*s0*s0),0.25) * cexp(I*k0*p.x[i])*cexp(-pow((p.x[i]-x0),2)/(2*s0*s0));		
}
