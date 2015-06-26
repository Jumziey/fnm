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
typedef struct par {
	double mgl;
	double i1;
	double i3;
}par;


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

//System equations
int
topLagrange(double time, const double v[], double d[], void *params)
{
	par *p = (par*)params;

	double mgl = p->mgl;
	double i1 = p->i1;
	double i3 = p->i3;
	
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
	d[5] = mgl/i1 * s + pow(v[3],2)*s*c - 
					i3/i1 * (v[4]+v[3]*c)*v[3]*s;
	
	return GSL_SUCCESS;
}

void
spin(double *v, double h, double t, double tEnd, double res, par *p, char* file, 
	int (*func)(double,const double*,double*,void*), double energy) 
{
	int status,i,steps;
	double k[2], verr[6], kDiff[2],tMeas, eDiff, c, s;
	gsl_odeiv2_system sys = {func, NULL, 6, p};
	gsl_odeiv2_step *stepMem = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 6);
	
	//Constants
	double i1 = p->i1;
	double i3 = p->i3;
	double mgl = p->mgl;
	s = sin(v[2]);c = cos(v[2]);
	k[0] = i1*v[3]*pow(s,2) + i3*(v[4]+v[3]*c)*c;
	k[1] = i3*(v[4]+v[3]*c);
	
	//Lets run simulation
	FILE *fp;
	fp = fopen(file,"w");
	tMeas = res;
	steps = ceil(abs((tEnd-t)/h));
	for(i=0;i<steps;i++){
		status = gsl_odeiv2_step_apply(stepMem, t, h, v, verr, NULL, NULL, &sys);
		if(status != GSL_SUCCESS){
			fprintf(stderr,"Problems!\n");
			break;
		}
		if(fabs(tMeas) > res) {
			c = cos(v[2]);s = sin(v[2]);
			kDiff[0] =k[0]- i1*v[3]*pow(s,2) - i3*(v[4]+v[3]*c)*c;
			kDiff[1] =k[1]- i3*(v[4]+v[3]*c);
			eDiff = energy- (i1/2 * (pow(v[5],2) + pow(v[3]*s,2)) + i3/2 * pow((v[4] + v[3]*c),2) + mgl*c);
			fprintf(fp,"%.6f %.6f %.6f \t %g %g %g\n", v[0], v[1], v[2], kDiff[0], kDiff[1], eDiff);
			tMeas = 0;
		}
		tMeas = tMeas+h;
		t = t+h;
	}
	fclose(fp);
	return;
}

int
main()
{
	int i;
	double abserr;
	size_t limit;
	par p;
	gsl_integration_workspace *w;
	gsl_function F;
	
	//We begin with integrating for the constants
	int funcs = 4;
	double consts[funcs];
	void *intFunc[] = {&mFunc, &mlFunc, &i1Func, &i3Func};
	
	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);
	
	for(i=0; i<funcs; i++) {
		F.function = intFunc[i];
		F.params = NULL;
		gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &consts[i] , &abserr);
	}
	//M, Ml, i1, i3
	p.mgl = consts[1]*9.82;
	p.i1 = consts[2];
	p.i3 = consts[3];

	//Here we run the simulation
	double v[] = {0.0, 0.0, M_PI/9, 0.0, 20*M_PI, 0.0};
	//Constants
	double i1 = p.i1;
	double i3 = p.i3;
	double mgl = p.mgl;
	double s = sin(v[2]);double c = cos(v[2]);
	double energy =  (i1/2 * (pow(v[5],2) + pow(v[3]*s,2)) + i3/2 * pow((v[4] + v[3]*c),2) + mgl*c);
	
	double h = 1e-5;
	double res = 2e-3;
	spin(v, h, 0, 4, res, &p,"./plots/data/forwardData", &topLagrange, energy);
	spin(v, -h, 4, 0, res, &p,"./plots/data/backwardData", &topLagrange, energy);

	return 0;
}
