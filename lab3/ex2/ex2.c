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
topLagrange(double t, const double v[], double d[], void *params)
{
	d[0] = v[3];
	d[1] = v[4];
	d[2] = v[5];
	
	d[3] = i3/i1*( v[3]*v[4]*(1/tan(v[2])) + v[4]*v[5]/sin(v[2]) )-2*v[3]*v[5];
	d[4] = v[3]*v[5]*sin(v[2]) + 2*v[3]*v[5]*pow(cos(v[2]),2)/sin(v[2]) - 
					i3/i1 * ( v[3]*v[5]*pow(cos(v[2]),2)/sin(v[2]) + v[4]*v[5]*(1/tan(v[2])) );
	d[5] = M*g*l/i1 * sin(v[2]) + pow(v[3],2)*sin(v[2])*cos(v[2]) - 
					i3/i1 * (v[4]+v[3]*cos(v[2]))*v[3]*sin(v[2]);
	
	return GSL_SUCCESS;
}

int
main()
{
	int funcs = 4;


	int i;
	double abserr, k[2], c[funcs];
	size_t limit;
	gsl_integration_workspace *w;
	gsl_function F;
	
	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);


	void *intFunc[] = {&mFunc, &mlFunc, &i1Func, &i3Func};
	for(i=0; i<funcs; i++) {
		F.function = intFunc[i];
		F.params = NULL;
		gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &c[i] , &abserr);
	}
	//M, l, i1, i3
	M = c[0];
	l = c[1] /= c[0];
	i1 = c[2];
	i3 = c[3];
	//+g!
	g = 9.81;
	printf("ans: %f\n", tan(2));
	
	//Now lets do some ode solving
	int status;
	double t = 0.0;
	double tEnd = 4.0;
	double h = 1e-5;
	double ti;
	gsl_odeiv2_system sys = {topLagrange, NULL, 6, NULL};
	gsl_odeiv2_driver *d = 
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h,epsabs, epsrel);
	
	double v[] = {0.0, 0.0, 20.0, 0.0, 10.0, 0.0};
	while(t<tEnd) {
		ti = t+h;
		status = gsl_odeiv2_driver_apply(d, &t, ti, v);
		if(status != GSL_SUCCESS){
			fprintf(stderr,"Problems!\n");
			break;
		} 
		k[0] = i1*v[3]*pow(sin(v[2]),2) + i3*(v[4]+v[3]*cos(v[2]))*cos(v[2]);
		k[1] = i3*(v[4]+v[3]*cos(v[2]));
		printf("%.6f %.6f %.6f \t%.6f %.6f\n", v[0], v[1], v[2], k[0], k[1]);
	}
	return 0;
}
