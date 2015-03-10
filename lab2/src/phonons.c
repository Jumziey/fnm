#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Need to declare
void frequencies(double A, double B, double m, double *q,
				double *omega, double *eps);

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
printVal(double *q, double *val) {
	int i;
	for(i=0; i<3; i++)
		printf("%g ", q[i]);
	for(i=0; i<3; i++)
		printf("%g ", val[i]);
	printf("\n");
}



void
freqEval(sp sub, double* q)
{
	double A,B, pdir[9], omega[3];
	
	//Calculating constants A and B from the problem specification
	A = sub.eps*( 156.0*(pow(sub.sigma,12)/pow(sub.rnn,14)) 
				- 84.0*(pow(sub.sigma,6)/pow(sub.rnn,8)) );
	B = 12.0* ( ( ( pow(sub.sigma,6)*(pow(sub.rnn,6)-pow(sub.sigma,6)) )/pow(sub.rnn,14) ) *sub.eps);
	
	frequencies(A, B, sub.m, q, omega, pdir);
	
	printVal(q,omega);
}

void
nFreqEval(sp sub, double* q1, double *q2, int n)
{
	int i,j;
	double diff[3], qCur[3];
	
	for(i=0; i<3; i++)
		diff[i] = (q2[i]-q1[i])/((double)n-1.0);
		
	for(i=0; i<n; i++) {
		for(j=0; j<3; j++)
			qCur[j] = q1[j]+diff[j]*i;
		freqEval(sub, qCur);
	}
}

void
volDepEval(sp sub, double *q)
{
	int i;
	//minus,normal, plus
	double omega[3][3], pdir[3][9], rnn[3], gamma[3], A[3], B[3];
	double h = 1E-20;// small compared to rnn
	
	//Initialize rnn values early for easy debug
	for(i=-1; i<2; i++)
		rnn[i+1] = sub.rnn + h*i;
	
	//We only have rnn dependency in A and B
	//(By extension only volume dependence here)
	for(i=0; i<3; i++){
		A[i] = sub.eps*( 156.0*(pow(sub.sigma,12)/pow(rnn[i],14)) 
					- 84.0*(pow(sub.sigma,6)/pow(rnn[i],8)) );
		B[i] = 12.0* ( ( ( pow(sub.sigma,6)*(pow(rnn[i],6)-pow(sub.sigma,6)) )/pow(rnn[i],14) ) *sub.eps);
	}
	
	for(i=0; i<3; i++) 
		frequencies(A[i], B[i], sub.m, q, omega[i], pdir[i]);
	
	//Array order is minus, normal, plus
	for(i=0; i<3; i++) {
		gamma[i] = (omega[2][i] - omega[0][i])/(pow(rnn[2],3)-pow(rnn[0],3));
		gamma[i] = fabs(gamma[i]*pow(rnn[1],3)/omega[1][i]);
	}
	printVal(q, gamma);
}
		
void
nValDepEval(sp sub, double *q1, double *q2, int n)
{
	int i,j;
	double diff[3], qCur[3];
	
	for(i=0; i<3; i++)
		diff[i] = (q2[i]-q1[i])/((double)n-1.0);
		
	for(i=0; i<n; i++) {
		for(j=0; j<3; j++)
			qCur[j] = q1[j]+diff[j]*i;
		volDepEval(sub, qCur);
	}
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
	double q1[3], q2[3];
	//Index ordering is coupled with subInd()
	const sp subs[4] = {Ne, Ar, Kr, Xe};
		
	if((subindex = subInd(argv[1]))<0)
		error("Unkown (for the program) substance\n");
	
	if((mode = operMode(argv[2]))<0)
		error("Mode is not implemented\n");
	
	int k;
	//Dealing with the command line arguments
	switch(mode) {
		case 0: //Omega
			switch(argc) {
				case 6:
					//We only give frequency for one point
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					freqEval(subs[subindex], q1);
					break;
				case 9:
					//Two points, frequency need to be evaluated in 11 points
					//Equally spaced
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nFreqEval(subs[subindex], q1, q2, 11);
					break;
				case 10:
					//Two points, frequency needs to evaluated in <npoints> points
					//Equally spaced
					npoints = atof(argv[9]);
					if(npoints<2)
						error("Need to evaluate two or more points if a range is given\n");
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nFreqEval(subs[subindex], q1, q2, npoints);
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		case 1: // gamma
			switch(argc) {
				case 6:
					//Calculating volume dep. of phonon frequency in one point
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					volDepEval(subs[subindex], q1);
					break;
				case 9:
					//Two points, need to be spaced 11
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nValDepEval(subs[subindex], q1, q2, 11);
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
					//printf("CRAZY!\n");
					break;
				case 5:
					//printf("blaoeuhtn\n");
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
			k=0;
			//printf("CongratZ! You've reach what the programmer thought to be an unreachable part of the program\n");		
	}
	
	return 0;	
}
