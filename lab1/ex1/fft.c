#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

//Damn c99!!!
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define SIZE 1024


double dt = 1.0/SIZE;

void
saveGrid(double* grid, int size,char* realFile, char* imFile) {
	int i;
	FILE* fp;	
	
	
	if(realFile) {
		fp = fopen(realFile, "w");
		for(i=0; i<size; i++)
			fprintf(fp,"%g\n",REAL(grid,i));
		fclose(fp);
	}
	if(imFile) {
		fp = fopen(imFile, "w");
		for(i=0; i<size; i++)
			fprintf(fp, "%g\n", IMAG(grid,i));
		fclose(fp);
	}
	
}

int main(){
	int i;
	double grid[2*SIZE];
	double alpha = 16.0*dt;
	double exponent = -dt*dt/(alpha*alpha);
	double c = 1.0/sqrt(alpha*alpha*M_PI);
	
	for(i=0; i<SIZE; i++) {
		REAL(grid,i) = c*exp(exponent*i*i);
		IMAG(grid,i) = 0;
	}
	fprintf(stderr, "choo!\n");
		
	saveGrid(grid,SIZE,"init", NULL);

	//gsl_fft_complex_radix2_forward(grid, 1, SIZE);
	gsl_fft_complex_radix2_dif_forward(grid,1,SIZE);
	//if( gsl_fft_complex_radix2_forward(grid, 1, SIZE) != GSL_SUCCESS) {
	//	fprintf(stderr, "Something went awry with the fourier transformation");
	//}
	saveGrid(grid, SIZE, "realTrans", "imTrans");
	}

