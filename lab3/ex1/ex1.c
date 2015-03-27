#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "top.h"

const double epsrel = 1E-12;
const double epsabs = 1E-8;

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
main()
{
	int funcs = 4;
	void *intFunc[] = {&mFunc, &mlFunc, &i1Func, &i3Func};

	int i;
	double m, ml, i1, i3, abserr, ans[funcs];
	size_t limit;
	gsl_integration_workspace *w;
	gsl_function F;
	
	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);
		
	for(i=0; i<funcs; i++) {
		F.function = intFunc[i];
		F.params = NULL;
		gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &ans[i] , &abserr);
	}
	printf("M: %.6f\n", ans[0]);
	printf("l: %.6f\n", ans[1]/ans[0]);
	printf("i1: %.6f\n", ans[2]);
	printf("i3: %.6f\n", ans[3]);
	
	
	return 0;
	
	
}
