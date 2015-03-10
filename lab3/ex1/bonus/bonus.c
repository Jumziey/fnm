#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "top.h"

const double epsrel = 1E-12;
const double epsabs = 1E-8;


int
main()
{
	int i;
	double n,max;
	
	n = 10000.0;
	max = 5.1;
	
	for(i=0; i<n; i++)
		printf("%g\n", top_y(max/n  * (double)i));
	return 0;	
}
