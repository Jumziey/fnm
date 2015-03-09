#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "frequencies.c"

//Substance Properties
typedef struct sp{
	double sigma;
	double eps;
	double rnn;
	double m;
}sp;

const sp Ne = {3.035E-10, 0.0721E-20, 3.1562E-10, 0.335092E-25};
const sp Ar = {3.709E-10, 0.236E-20, 3.7477E-10, 0.66335E-25};
const sp Kr = {3.966E-10, 0.325E-20, 3.9922E-10, 1.3915E-25};
const sp Xe = {4.318E-10, 0.458E-20, 4.3346E-10, 2.18017E-25};
;

void
error(char* err) {
	fprintf(stderr,"%s", err);
	exit(2);
}

void
printOmegaVal(double *q, double *omega) {
	int i;
	for(i=0; i<3; i++)
		printf("%g ", q[i]);
	for(i=0; i<3; i++)
		printf("%g ", omega[i]);
	printf("\n");
}



void
pointEval(sp sub, double* q)
{
	double A,B, *pdir, *omega;
	
	//Initialize to zero for frequencies()
	pdir = calloc(3*3,sizeof(double));
	omega = calloc(3,sizeof(double));
	
	//Calculating constants A and B from the problem specification
	A = sub.eps*( 156.0*(pow(sub.sigma,12)/pow(sub.rnn,14)) 
				- 84.0*(pow(sub.sigma,6)/pow(sub.rnn,8)) );
	B = 12.0* ( ( ( pow(sub.sigma,6)*(pow(sub.rnn,6)-pow(sub.sigma,6)) )/pow(sub.rnn,14) ) *sub.eps);
	
	frequencies(A, B, sub.m, q, omega, pdir);
	
	printOmegaVal(q,omega);
	free(pdir);
	free(omega);
}

void
nPointsEval(sp sub, double* q1, double *q2, int n)
{
	int i,j;
	double* diff, *qCur;
	diff = malloc(3*sizeof(double));
	qCur = malloc(3*sizeof(double));
	
	for(i=0; i<3; i++)
		diff[i] = (q2[i]-q1[i])/((double)n-1.0);
		
	for(i=0; i<n; i++) {
		for(j=0; j<3; j++)
			qCur[j] = q1[j]+diff[j]*i;
		pointEval(sub, qCur);
	}
	free(diff);
	free(qCur);
}
		



int
subInd(char* substance)
{
	if(!strncmp(substance, "Ne", 2))
		return 0;
	if(!strncmp(substance, "Ar", 2))
		return 1;
	if(!strncmp(substance, "Kr", 2))
		return 2;
	if(!strncmp(substance, "Xe", 2))
		return 3;
		
	return -1;
}

int
operMode(char* oper)
{
	if(!strncmp(oper, "omega",5))
		return 0;
	if(!strncmp(oper, "gamma",5))
		return 1;
	if(!strncmp(oper, "cv",2))
		return 2;
	
	return -1;
}
	 
int
main(int argc, char** argv)
{
	int mode, subindex, npoints;
	double *q1, *q2;
	//Index ordering is coupled with subInd()
	const sp subs[4] = {Ne, Ar, Kr, Xe};

	
		
	if((subindex = subInd(argv[1]))<0)
		error("Unkown (for the program) substance\n");
	
	if((mode = operMode(argv[2]))<0)
		error("Mode is not implemented\n");
	
	
	//Dealing with the command line arguments
	switch(mode) {
		case 0: //Omega
			switch(argc) {
				case 6:
					q1 = malloc(sizeof(double)*3);
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					//printf("only one point\n");
					pointEval(subs[subindex], q1);
					free(q1);
					break;
				case 9:
					q1 = malloc(sizeof(double)*3);
					q2 = malloc(sizeof(double)*3);
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nPointsEval(subs[subindex], q1, q2, 11);
					//Two points, need to be spaced 11
					free(q1);free(q2);
					break;
				case 10:
					npoints = atof(argv[9]);
					if(npoints<2)
						error("Need to evaluate two or more points if a range is given\n");
					q1 = malloc(sizeof(double)*3);
					q2 = malloc(sizeof(double)*3);
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nPointsEval(subs[subindex], q1, q2, npoints);
					free(q1);free(q2);
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		case 1: // gamma
			switch(argc) {
				case 6:
					printf("CRAZY!\n");
					break;
				case 9:
					printf("blaoeuhtn\n");
					//Two points, need to be spaced 11
					break;
				case 10:
					//Two points and npoints
					if(atoi(argv[9])<2)
						error("Need to evaluate 2 or more points\n"); 
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		case 2: // cv
			switch(argc) {
				case 4:
					printf("CRAZY!\n");
					break;
				case 5:
					printf("blaoeuhtn\n");
					//Two points, need to be spaced 11
					break;
				case 6:
					//Two points and npoints
					if(atoi(argv[5])<2)
						error("Need to evaluate 2 or more points\n"); 
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		default:
			printf("CongratZ! You've reach an unreachable part of the program\n");
	}
	
	
	


	return 0;	
}
