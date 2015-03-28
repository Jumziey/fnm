/*
This is a collection of debugging functions, created for this particular program
Should be usable for a collection of GSL programs
*/

void
print_matrix(gsl_matrix *m, int r, int c)
{
	int i,j;
	
	fprintf(stderr,"####################\n");
	for(i=0; i<c; i++) {
		for(j=0; j<c; j++) {
			fprintf(stderr,"%g \t", gsl_matrix_get(m, i, j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"####################\n");
}

void
printEig(double *omega, double *eps) {
	int i,j;
	fprintf(stderr,"qqqqqqqqqqqqqqqqqqqqqqqqqq\n");
	for(i=0; i<3; i++) {
		fprintf(stderr,"Eigenvalue: %g",omega[i]);
		fprintf(stderr,"\t Eigenvector:");
		for(j=0;j<3;j++)
			fprintf(stderr,"\t%g", eps[i*3+j]);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"qqqqqqqqqqqqqqqqqqqqqqqqqq\n");
}

