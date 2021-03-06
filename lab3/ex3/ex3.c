#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "top.h"

const double epsrel = 1E-14;
const double epsabs = 1E-10;

//System Constants
typedef struct par {
	double mgl;
	double i1;
	double i3;
}par;

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
topHam(double time, const double v[], double d[], void *params)
{
	par *p = (par*)params;

	double mgl = p->mgl;
	double i1 = p->i1;
	double i3 = p->i3;

	const double c = cos(v[2]);
	const double s = sin(v[2]);
	
	d[0] = (2.0*v[3] - 2.0*v[4]*c)/(2.0*i1*pow(s,2));
	d[1] = v[4]/i3 - (c*(v[3] - v[4]*c))/(i1*pow(s,2));
	d[2] = v[5]/i1;
	
	d[3] = 0;
	d[4] = 0;
	d[5] = mgl*s - (v[4]*(v[3] - v[4]*c))/(i1*s) + (c*pow((v[3] - v[4]*c),2))/(i1*pow(s,3));
	
	
	return GSL_SUCCESS;
}


int
topJac(double time, const double v[], double *dfdv, double *dfdt, void *params)
{
	int i;

	par *p = (par*)params;

	double mgl = p->mgl;
	double i1 = p->i1;
	double i3 = p->i3;
	double v0=v[0],v1=v[1],v2=v[2],v3=v[3],v4=v[4],v5=v[5];
	
	gsl_matrix_view dfdv_mat;
	gsl_matrix *m;
	
	const double ct = cos(v[2]);
 	const double st = sin(v[2]);
 	
	dfdv_mat = gsl_matrix_view_array(dfdv, 6, 6);
	m = &dfdv_mat.matrix; 
 	
	gsl_matrix_set(m,0,0,0);
	gsl_matrix_set(m,0,1,0);
	gsl_matrix_set(m,0,2,v4/(i1*st) - ct*(2.0*v3 - 2.0*v4*ct)/(i1*pow(st,3)));
	gsl_matrix_set(m,0,3,1/(i1*pow(st,2)));
	gsl_matrix_set(m,0,4,-ct/(i1*pow(st,2)));
	gsl_matrix_set(m,0,5,0);
	
	gsl_matrix_set(m,1,0,0);
	gsl_matrix_set(m,1,1,0);
	gsl_matrix_set(m,1,2,(v3 - v4*ct)/(i1*st) + (2.0*pow(ct,2)*(v3 - v4*ct))/(i1*pow(st,3)) - (v4*ct)/(i1*st));
	gsl_matrix_set(m,1,3,-ct/(i1*pow(st,2)));
	gsl_matrix_set(m,1,4,1/i3 + pow(ct,2)/(i1*pow(st,2)));
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
	gsl_matrix_set(m,5,2,mgl*ct - v4*st/i1 - pow((v3 - v4*ct),2)/(i1*pow(st,2)) - (3.0*pow(ct,2)*pow((v3 - v4*ct),2))/(i1*pow(st,4)) + (3.0*v4*ct*(v3 - v4*ct))/(i1*pow(st,2)));
	
	gsl_matrix_set(m,5,3,(ct*(2.0*v3 - 2.0*v4*ct))/(i1*pow(st,3)) - v4/(i1*st));
	
	gsl_matrix_set(m,5,4,(v4*ct)/(i1*st) - (2.0*pow(ct,2)*(v3 - v4*ct))/(i1*pow(st,3)) - (v3 - v4*ct)/(i1*st));
	gsl_matrix_set(m,5,5,0);
	
	for(i=0;i<6;i++)
		dfdt[i] = 0;
		
	return GSL_SUCCESS;
}

void
spin(double *v, double h, double t, double tEnd, double res, par *p, char* file, 
	int (*func)(double,const double*,double*,void*),
	int (*jac) (double,const double*,double*,double*,void*), double energy) 
{
	int status,i,steps;
	double  verr[6],tMeas, eDiff;
	gsl_odeiv2_system sys = {func, jac, 6, p};
	gsl_odeiv2_step *stepMem = gsl_odeiv2_step_alloc(gsl_odeiv2_step_bsimp, 6);
	
	//Constants
	double i1 = p->i1;
	double i3 = p->i3;
	double mgl = p->mgl;
	double c = cos(v[2]), s = sin(v[2]);

		
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
			c = cos(v[2]); s = sin(v[2]);
			eDiff = energy - (pow(v[5],2)/(2*i1) + pow((v[3]-v[4]*c),2)/(2*i1*pow(s,2)) + pow(v[4],2)/(2*i3) + mgl*c);
			fprintf(fp,"%.8f %.8f %.8f %e %e %e \n", v[0], v[1], v[2], v[3], v[4], eDiff);
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
	//Mgl, i1, i3
	p.mgl = consts[1]*9.82;
	p.i1 = consts[2];
	p.i3 = consts[3];

	//Here we run the simulation
	//double v[] = {0.0, 0.0, M_PI/9, p.i3*20*M_PI*cos(M_PI/9), p.i3*20*M_PI, 0.0};
	double v[] = {0.0, 0.0, M_PI/9, p.i3*cos(M_PI/9)*20*M_PI, p.i3*20*M_PI, 0.0};
	double c = cos(v[2]), s = sin(v[2]);
	double energy = pow(v[5],2)/(2*p.i1) + pow((v[3]-v[4]*c),2)/(2*p.i1*pow(s,2)) + pow(v[4],2)/(2*p.i3) + p.mgl*c;
	double h = 1e-5;
	double res = 2e-3;
	spin(v, h, 0, 4, res, &p,"./plots/data/forwardData", &topHam, &topJac, energy);
	spin(v, -h, 4, 0, res, &p,"./plots/data/backwardData", &topHam, &topJac, energy);

	return 0;
}

