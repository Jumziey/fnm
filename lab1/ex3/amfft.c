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
double
pwrMean(double* grid) {
	int i=0;
	double sum = 0.0;
	for(i=0; i<SIZE; i++)
		sum += REAL(grid,i)*REAL(grid,i)+IMAG(grid,i)*IMAG(grid,i);
	return(sum/(double)SIZE);
}

int main(int argc, char* argv[]){
	int i;
	double  grid[2*SIZE],grid2[2*SIZE];
	FILE* fp;
	char line[60];
	
	fp = fopen(argv[1], "r");
	for(i=0; i<SIZE; i++) {
		fgets(line, 60, fp);
		REAL(grid,i) = atof(line);
		IMAG(grid,i) = 0;
	}
	fclose(fp);
	
	gsl_fft_complex_radix2_forward(grid, 1, SIZE);
	saveGrid(grid, SIZE, "realTrans", "imagTrans");
	
	//Apply filter
	double fc = 1024.0;
	double bw = 256.0;
	double exponent;
	for(i=0; i<SIZE/2-1; i++) {
		exponent = -0.5*(((double)i - fc)/bw)*((i - fc)/bw);
		REAL(grid,i) *= expf(exponent);
		IMAG(grid,i) *= expf(exponent);
		REAL(grid,SIZE-i) *= expf(exponent);
		IMAG(grid,SIZE-i) *= expf(exponent);
	}
	REAL(grid,SIZE/2) = 0; //We should have zero at frequency 0
	IMAG(grid,SIZE/2) = 0;
	saveGrid(grid,SIZE, "realFilteredTrans", "imagFilteredTrans");
	
	//Taking back the filtered signal to time space
	gsl_fft_complex_radix2_backward(grid,1,SIZE);
	for(i=0; i<SIZE; i++)
		REAL(grid,i) *= 1.0/SIZE; // Normalizing
	saveGrid(grid,SIZE, "filteredSignal", NULL);
	
	//Decoding the bit message
	double bit[SIZE/64];
	for(i=0; i<SIZE; i++) {
		bit[i/64] += REAL(grid,i)*REAL(grid,i)+IMAG(grid,i)*IMAG(grid,i);
	}
	double mean = pwrMean(grid)*64;
	char message[SIZE/64/8+1];
	
	for(i=0; i<SIZE/64; i++) {
		printf("bit[i] = %g, mean = %g\n",bit[i]/64,mean);
		message[i/8] = message[i/8]*2 + (bit[i] > mean ? 1:0);
	}
	
	printf("%s\n",message);
	
	return 0;
}

