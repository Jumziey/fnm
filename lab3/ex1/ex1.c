#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "top.h"

const double epsrel = 1E-12;
const double epsabs = 1E-8;

double
top(double x, void* params)
{
	return M_PI*top_y(x)*top_y(x)*top_rho(x);
}

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
	
int
main()
{
	double M, l, I1, I3, abserr;
	size_t limit;
	gsl_integration_workspace *w;
	gsl_function F;
	F.function = &top;
	F.params = NULL;

	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);
	gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &M, &abserr);
	
	printf("M: %g\n", M);
	
	return 0;
	
	
}
