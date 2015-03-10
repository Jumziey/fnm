#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "top.h"

const double epsrel = 1E-12;
const double epsabs = 1E-8;

double
cyl(double x, void* params)
{
	return top_y(x)*top_y(x)*top_rho(x);
}
int
main()
{
	double result, abserr;
	size_t limit;
	gsl_integration_workspace *w;
	gsl_function F;
	F.function = &cyl;
	F.params = NULL;

	limit = 1000;
	w = gsl_integration_workspace_alloc(limit);
	gsl_integration_qags(&F, 0, 5, epsabs, epsrel, limit, w, &result, &abserr);
	
	printf("result: %g\n", result);
	
	return 0;
	
	
}
