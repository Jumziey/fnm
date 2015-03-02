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
	int u0 = 1;
	int u1 = 3;
	double Fc = (1.0/(128.0*dt))*2;
	
	int u = u1; 
	for(i=0; i<SIZE; i++) {
		if(i%(SIZE/4) == 0) {
			if(u == u0)
				u = u1;
			else
				u = u0;
		}
		REAL(grid,i) = u * sin(2*M_PI * Fc * i * dt);
		IMAG(grid,i) = 0;
	}
	
	//Saving before and after Fourier transform
	saveGrid(grid,SIZE,"init", NULL);
	gsl_fft_complex_radix2_dif_forward(grid,1,SIZE);
	saveGrid(grid, SIZE, "realTrans", "imTrans");
	}

