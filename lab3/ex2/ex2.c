#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "top.h"

const double epsrel = 1E-12;
const double epsabs = 1E-8;

//System Constants
//Yes i could use params... But LAZY! :D
double M,l,g,i1,i3;


//For reverse engineering
void
saveData(double (*func)(double), char* file, double s, double e, int n)
{
	int i;
	double diff;
	FILE* fp;
	
	fp = fopen(file, "w");
	diff = (e-s)/(double)n;
	for(i=0; i<n; i++)
		fprintf(fp,"%g %g\n", s+diff*i, func(s+diff*i));
	fclose(fp);
}

//Integral expressions
double
mFunc(double x, void* params)
{
	return M_PI*pow(top_y(x
	),2)*top_rho(x);
}

double
mlFunc(double x, void* params)
{
	return  M_PI*pow(top_y(x),2)*top_rho(x)*x;
}

double
i1Func(double x, void* params)
{
	return M_PI * top_rho(x)*pow(top_y(x),2)*pow(x,2)
				+ M_PI/4 * top_rho(x)*pow(top_y(x),4);
}

double
i3Func(double x, void* params)
{
	return M_PI/2 * top_rho(x)*pow(top_y(x),4);
}


int
topLagrange(double time, const double v[], double d[], void *params)
{
	const double c = cos(v[2]);
	const double s = sin(v[2]);
	const double t = tan(v[2]);
	const double ct = 1.0/tan(v[2]);
	
	d[0] = v[3];
	d[1] = v[4];
	d[2] = v[5];
	
	d[3] = i3/i1*( v[3]*v[5]*ct + v[4]*v[5]/s )-2*v[3]*v[5]*ct;
	
	
	d[4] = v[3]*v[5]*s + 2*v[3]*v[5]*ct*c - 
					i3/i1 * ( v[3]*v[5]*pow(c,2)/s + v[4]*v[5]*ct );
					
	d[5] = M*g*l/i1 * s + pow(v[3],2)*s*c - 
					i3/i1 * (v[4]+v[3]*c)*v[3]*s;
	
	return GSL_SUCCESS;
}

int
main()
{
	int funcs = 4;


	int i;
	double abserr, consts[funcs];
	size_t limit;
	gsl_integration_workspace *w;
	gsl_function F;
	
	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);


	void *intFunc[] = {&mFunc, &mlFunc, &i1Func, &i3Func};
	for(i=0; i<funcs; i++) {
		F.function = intFunc[i];
		F.params = NULL;
		gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &consts[i] , &abserr);
	}
	//M, l, i1, i3
	M = consts[0];
	l = consts[1] /= consts[0];
	i1 = consts[2];
	i3 = consts[3];
	g = 9.81;
	
	//Now lets do some ode solving
	int status;
	double k[2], kDiff[2], energy, eDiff, c, s;
	double t = 0.0;
	double tEnd = 4.0;
	double res = 1e-4;
	double h = 1e-4;
	double ti;
	gsl_odeiv2_system sys = {topLagrange, NULL, 6, NULL};
	gsl_odeiv2_driver *d = 
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h,epsabs, epsrel);
	
	//initial values
	double v[] = {0.0, 0.0, 20.0, 0.0, 10.0, 0.0};
	k[0] = i1*v[3]*pow(sin(v[2]),2) + i3*(v[4]+v[3]*cos(v[2]))*cos(v[2]);
	k[1] = i3*(v[4]+v[3]*cos(v[2]));
	printf("%g %g\n", k[0], k[1]);
	

	s = sin(v[2]);c = cos(v[2]);
	energy =  (i1/2 * (v[5]*v[5] + v[3]*v[3]*s) + i3/2 * pow((v[4] + v[3]*c),2) + M*g*l*c);
	printf("energy: %g\n", energy);
	while(t<tEnd) {
		ti = t+res;
		status = gsl_odeiv2_driver_apply(d, &t, ti, v);
		if(status != GSL_SUCCESS){
			fprintf(stderr,"Problems!\n");
			break;
		}
		c = cos(v[2]);s = sin(v[2]);
		kDiff[0] =k[0]- i1*v[3]*pow(sin(v[2]),2) - i3*(v[4]+v[3]*c)*c;
		kDiff[1] =k[1]- i3*(v[4]+v[3]*cos(v[2]));
		eDiff = energy- (i1/2 * (v[5]*v[5] + v[3]*v[3]*s) + i3/2 * pow((v[4] + v[3]*c),2) + M*g*l*c);
		printf("%.6f %.6f %.6f \t %g %g %g\n", v[0], v[1], v[2], kDiff[0], kDiff[1], eDiff);
	}
	return 0;
}
