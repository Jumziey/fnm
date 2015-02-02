#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_matrix.h>

//Damn c99!!!
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#define SIZE 1024

double dt = 1/SIZE;

void
printArray(double* grid, int size) {
	int i;
	
	for(i=0; i<size; i++)
		printf("%f\n",grid[i]);
}

int main(){
	int i;
	double grid[SIZE];
	long double alpha = 16*dt;
	printf("alpha: %.4f\n",alpha);
	long double exponent;
	long double c = 1/(sqrt(M_PI*alpha*alpha));
	printf("const: %g\n",c);
	return 0;
	for(i=0; i<SIZE; i++) {
		exponent = i*i*dt*dt;
		
		printf("%g\n",exp(-exponent));
	}
}

