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

#define SIZE 8192


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
	FILE* fp;
	char line[60];
	
	fp = fopen("am26.dat", "r");
	
	for(i=0; i<SIZE; i++) {
		fgets(line, 60, fp);
		REAL(grid,i) = atof(line);
		//fgets(line, 60, fp);
		IMAG(grid,i) = 0;
	}
		
	saveGrid(grid,SIZE,"init", NULL);
	gsl_fft_complex_radix2_forward(grid, 1, SIZE);
	saveGrid(grid, SIZE, "realTrans", NULL);
	//Now to filter the data in grid
	double fc = 1152;
	double bw = 1280-895.6;
	//Positive part 0 - max
	for(i=0; i<SIZE/2; i++) 
		REAL(grid,i) = expf( (-1/2) * powf(((1.000122085215480*i - fc)/bw),2) );
	//Negative part 0 - min
	for(i=SIZE/2; i<SIZE; i++)
		REAL(grid,i) = expf( (-1/2) * powf(((1.000122085215480*(i-SIZE/2+1) - fc)/bw),2) );
	
	gsl_fft_complex_radix2_backward(grid, 1, SIZE);
	
	saveGrid(grid, SIZE, "filteredSignal", NULL);
	}

