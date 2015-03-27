#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

const double epsrel = 1E-12;
const double epsabs = 1E-8;

double
s(double r)
{
	return 0.5*exp(-r)*r;
}

double
phiAn(double r)
{
	return 1-0.5*(r+2.0)*exp(-r);
}

void
iterate(double h, double phi, char *file, char *table)
{
	int i;
	double r,rEnd,phim,phip,err,j;
	FILE *fp1, *fp2;
	
	j=1;
	h = 1e-1;
	rEnd = 20;
	fp1 = fopen(file,"w");
	fp2 = fopen(table,"w");
	phim = 0;
	r = h;
	fprintf(fp1,"%g %g 0\n", r-h,phim);
	fprintf(fp1,"%g %g %g\n", r,phi, phiAn(r)-phi);
	for(;r<rEnd;r+=h) {
		phip = 2.0*phi - (h*h/12.0)*(s(r+h)+10.0*s(r)+s(r-h))-phim;
		phim = phi;
		phi = phip;
		err = phiAn(r+h)-phi;
		fprintf(fp1,"%g %g %g\n", r+h,phi, err);
		if(fabs((r+h)-(2*j))<1e-10) {
			fprintf(fp2,"%g & %g & %g \\\\ \\hline \n",r+h, phi, err);
			j++;
		}
	}
	fclose(fp1);
	fclose(fp2);
	return;
}

int
main()
{
	double h, phi, abserr;
	
	h = 1e-1;
	//analytical
	iterate(h, phiAn(h), "./plots/data/phi1", "./plots/data/phi1table");
	
	//Numerical Integration
	size_t limit = 1000;
	gsl_function F;
	gsl_integration_workspace *w;
	
	w = gsl_integration_workspace_alloc(limit);
	F.function = &s;
	F.params = NULL;
	gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w, &phi, &abserr);
	iterate(h, phi*h, "./plots/data/phi2", "./plots/data/phi2table");
	
	//On the 95% limit
	iterate(h,phiAn(h)*0.95, "./plots/data/phi3", "./plots/data/phi3table");
	
	return 0;
}
















	
