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
topHam(double time, const double v[], double d[], void *params)
{
	const double c = cos(v[2]);
	const double s = sin(v[2]);
	
	d[0] = (2*v[3] - 2*v[4]*c)/(2*i1*s*s);
	d[1] = v[4]/i3 - (c*(v[3] - v[4]*c))/(i1*s*s);
	d[2] = v[5]/i1;
	
	d[3] = 0;
	d[4] = 0;
	d[5] = M*l*g*s - (v[4]*(v[3] - v[4]*c))/(i1*s) + (c*pow((v[3] - v[4]*c),2))/(i1*s*s*s);
	
	return GSL_SUCCESS;
}

int
topJac(double time, const double v[], double *dfdv, double *dfdt, void *params)
{
	int i;
	gsl_matrix_view dfdv_mat;
	gsl_matrix *m;
	
	const double c = cos(v[2]);
 	const double s = sin(v[2]);
 	
	dfdv_mat = gsl_matrix_view_array(dfdv, 6, 6);
	m = &dfdv_mat.matrix; 
 	
	gsl_matrix_set(m,0,0,0);
	gsl_matrix_set(m,0,1,0);
	gsl_matrix_set(m,0,2,v[4]/(i1*s) - (c*(2*v[3] - 2*v[4]*c))/(i1*pow(s,3)));
	gsl_matrix_set(m,0,3,1/(i1*s*s));
	gsl_matrix_set(m,0,4,-c/(i1*s*s));
	gsl_matrix_set(m,0,5,0);
	gsl_matrix_set(m,1,0,0);
	gsl_matrix_set(m,1,1,0);
	gsl_matrix_set(m,1,2,(v[3] - v[4]*c)/(i1*s) - (v[4]*c)/(i1*s) + (2*c*c*(v[3] - v[4]*c))/(i1*s*s*s));
	gsl_matrix_set(m,1,3,-c/(i1*s*s));
	gsl_matrix_set(m,1,4,1/i3 + c*s/(i1*s*s));
	gsl_matrix_set(m,1,5,0);
	gsl_matrix_set(m,2,0,0);
	gsl_matrix_set(m,2,1,0);
	gsl_matrix_set(m,2,2,0);
	gsl_matrix_set(m,2,3,0);
	gsl_matrix_set(m,2,4,0);
	gsl_matrix_set(m,2,5,1/i1);
	gsl_matrix_set(m,3,0,0);
	gsl_matrix_set(m,3,1,0);
	gsl_matrix_set(m,3,2,0);
	gsl_matrix_set(m,3,3,0);
	gsl_matrix_set(m,3,4,0);
	gsl_matrix_set(m,3,5,0);
	gsl_matrix_set(m,4,0,0);
	gsl_matrix_set(m,4,1,0);
	gsl_matrix_set(m,4,2,0);
	gsl_matrix_set(m,4,3,0);
	gsl_matrix_set(m,4,4,0);
	gsl_matrix_set(m,4,5,0);
	gsl_matrix_set(m,5,0,0);
	gsl_matrix_set(m,5,1,0);
	gsl_matrix_set(m,5,2,M*l*g*c - v[4]*v[4]/i1 - pow((v[3] - v[4]*c),2)/(i1*s*s) - (3*c*c*pow((v[3] - v[4]*c),2))/(i1*pow(s,4)) + (3*v[4]*c*(v[3] - v[4]*c))/(i1*s*s));
	gsl_matrix_set(m,5,3,(c*(2*v[3] - 2*v[4]*c))/(i1*s*s*s) - v[4]/(i1*s));
	gsl_matrix_set(m,5,4,(v[4]*c)/(i1*s) - (v[3] - v[4]*c)/(i1*s) - (2*c*c*(v[3] - v[4]*c))/(i1*s*s*s));
	gsl_matrix_set(m,5,5,0);
	
	for(i=0;i<6;i++)
		dfdt[i] = 0;
		
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
	double h = 2e-3;
	double ti;
	gsl_odeiv2_system sys = {topHam, topJac, 6, NULL};
	gsl_odeiv2_driver *d = 
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp, h,epsabs, epsrel);
	
	//initial values
	double v[] = {0.0, 0.0, 20.0, 0.0, 10.0, 0.0};
	k[0] = i1*v[3]*pow(sin(v[2]),2) + i3*(v[4]+v[3]*cos(v[2]))*cos(v[2]);
	k[1] = i3*(v[4]+v[3]*cos(v[2]));

	s = sin(v[2]);c = cos(v[2]);
	energy =  (i1/2 * (pow(v[5],2) + pow(v[3]*s,2)) + i3/2 * pow((v[4] + v[3]*c),2) + M*g*l*c);
	
	//Running simulation forwards
	FILE *fp;
	fp = fopen("./plots/data/forwardData","w");
	while(t<tEnd) {
		ti = t+h;
		status = gsl_odeiv2_driver_apply(d, &t, ti, v);
		if(status != GSL_SUCCESS){
			fprintf(stderr,"Problems!\n");
			break;
		}
		c = cos(v[2]);s = sin(v[2]);
		kDiff[0] =k[0]- i1*v[3]*pow(sin(v[2]),2) - i3*(v[4]+v[3]*c)*c;
		kDiff[1] =k[1]- i3*(v[4]+v[3]*cos(v[2]));
		eDiff = energy- (i1/2 * (pow(v[5],2) + pow(v[3]*s,2)) + i3/2 * pow((v[4] + v[3]*c),2) + M*g*l*c);
		fprintf(fp,"%.6f %.6f %.6f \t %g %g %g\n", v[0], v[1], v[2], kDiff[0], kDiff[1], eDiff);
	}
	fclose(fp);
	
	//Running simulation backwards
	gsl_odeiv2_driver_free(d);
	d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, -h,epsabs, epsrel);
	tEnd = 0;
	fp = fopen("./plots/data/backwardData","w");
	while(t>0) {
		ti = t-h;
		status = gsl_odeiv2_driver_apply(d, &t, ti, v);
		if(status != GSL_SUCCESS){
			fprintf(stderr,"Problems!\n");
			break;
		}
		c = cos(v[2]);s = sin(v[2]);
		kDiff[0] =k[0]- i1*v[3]*pow(sin(v[2]),2) - i3*(v[4]+v[3]*c)*c;
		kDiff[1] =k[1]- i3*(v[4]+v[3]*cos(v[2]));
		eDiff = energy- (i1/2 * (pow(v[5],2) + pow(v[3]*s,2)) + i3/2 * pow((v[4] + v[3]*c),2) + M*g*l*c);
		fprintf(fp,"%.6f %.6f %.6f \t %g %g %g\n", v[0], v[1], v[2], kDiff[0], kDiff[1], eDiff);
	}
	fclose(fp);
	return 0;
}
