#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

//I do declare
void frequencies(double A, double B, double m, double *q,
				double *omega, double *eps);

//Substance Properties
typedef struct sp{
	double sigma;
	double eps;
	double rnn;
	double m;
}sp;

//Defining constants
const double kb = 1.3806488E-23;
const double hbar = 1.0545717E-34;

//Substances
const sp Ne = {3.035E-10, 0.0721E-20, 3.1562E-10, 0.335092E-25};
const sp Ar = {3.709E-10, 0.236E-20, 3.7477E-10, 0.66335E-25};
const sp Kr = {3.966E-10, 0.325E-20, 3.9922E-10, 1.3915E-25};
const sp Xe = {4.318E-10, 0.458E-20, 4.3346E-10, 2.18017E-25};

//What!? You wanted me to read the file!? Screw that! 
//I'm waay to lazy for that ;)
double qw[] = {	 1.000,0.400,0.000,12.00
				,1.000,0.200,0.200,12.00
				,1.000,0.200,0.000,12.00
				,1.000,0.000,0.000, 3.00
				,0.900,0.500,0.100,24.00
				,0.900,0.300,0.300,12.00
				,0.900,0.300,0.100,48.00
				,0.900,0.100,0.100,24.00
				,0.800,0.600,0.000,24.00
				,0.800,0.400,0.200,48.00
				,0.800,0.400,0.000,24.00
				,0.800,0.200,0.200,24.00
				,0.800,0.200,0.000,24.00
				,0.800,0.000,0.000, 6.00
				,0.700,0.700,0.100,12.00
				,0.700,0.500,0.300,24.00
				,0.700,0.500,0.100,48.00
				,0.700,0.300,0.300,24.00
				,0.700,0.300,0.100,48.00
				,0.700,0.100,0.100,24.00
				,0.600,0.600,0.200,24.00
				,0.600,0.600,0.000,12.00
				,0.600,0.400,0.400,24.00
				,0.600,0.400,0.200,48.00
				,0.600,0.400,0.000,24.00
				,0.600,0.200,0.200,24.00
				,0.600,0.200,0.000,24.00
				,0.600,0.000,0.000, 6.00
				,0.500,0.500,0.500, 4.00
				,0.500,0.500,0.300,24.00
				,0.500,0.500,0.100,24.00
				,0.500,0.300,0.300,24.00
				,0.500,0.300,0.100,48.00
				,0.500,0.100,0.100,24.00
				,0.400,0.400,0.400, 8.00
				,0.400,0.400,0.200,24.00
				,0.400,0.400,0.000,12.00
				,0.400,0.200,0.200,24.00
				,0.400,0.200,0.000,24.00
				,0.400,0.000,0.000, 6.00
				,0.300,0.300,0.300, 8.00
				,0.300,0.300,0.100,24.00
				,0.300,0.100,0.100,24.00
				,0.200,0.200,0.200, 8.00
				,0.200,0.200,0.000,12.00
				,0.200,0.000,0.000, 6.00
				,0.100,0.100,0.100, 8.00
				,0.025,0.020,0.015, 1.00};

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

double*
freqEval(sp sub, double* q)
{
	double A,B, pdir[9], *omega;
	
	//We will return omega
	omega = calloc(3,sizeof(double));
	
	//Calculating constants A and B from the problem specification
	A = sub.eps*( 156.0*(pow(sub.sigma,12)/pow(sub.rnn,14)) 
				- 84.0*(pow(sub.sigma,6)/pow(sub.rnn,8)) );
	B = 12.0* ( ( ( pow(sub.sigma,6)*(pow(sub.rnn,6)-pow(sub.sigma,6)) )/pow(sub.rnn,14) ) *sub.eps);
	
	frequencies(A, B, sub.m, q, omega, pdir);
	return omega;
}

double*
volDepEval(sp sub, double *q)
{
	int i;
	//minus,normal, plus
	double omega[3][3], pdir[3][9], rnn[3], A[3], B[3], *gamma;
	double h = 1E-20;// small compared to rnn
		
	gamma = calloc(3, sizeof(double));
	
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
		//Can't divide by zero, that would be silly from all kinds of perspectives.
		if(omega[1][i] == 0)
			return NULL;
		gamma[i] = (omega[2][i] - omega[0][i])/(pow(rnn[2],3)-pow(rnn[0],3));
		gamma[i] = fabs(gamma[i]*pow(rnn[1],3)/omega[1][i]);
	}
	return gamma;
}

void
nEval(sp sub, double *q1, double *q2, int n, double* (*evalFunc)(sp, double*))
{
	int i,j;
	double diff[3], qCur[3], *ret;
	
	for(i=0; i<3; i++)
		diff[i] = (q2[i]-q1[i])/((double)n-1.0);
		
	for(i=0; i<n; i++) {
		for(j=0; j<3; j++)
			qCur[j] = q1[j]+diff[j]*i;
		ret = evalFunc(sub, qCur);
		if(ret != NULL)
			printVal(qCur, ret);
	}
}

double
cvEval(sp sub, double T)
{
	int i,j;
	double pdir[9], *omega, Cv, H;
	
	Cv = 0;
	//We have 48 rows with 3 q values (x,y,z) and a lastly weight value W
	for(i=0; i<48; i++) {
		//First we calculate the frequency
		omega = freqEval(sub, &qw[i*4]); //Pointer Arithmetic! Gotta love :D
		for(j=0; j<3; j++) {
			H = hbar*omega[j]/(kb*T);
			Cv += qw[i*4 + 3] * kb * H*H * exp(H)/pow((exp(H)-1),2);
		}
	}
	Cv *= 4*pow(M_PI,3)/(pow(M_PI,3)*8000*pow(sub.rnn/sqrt(2),3));
	return Cv;
}
	
void
nCvEval(sp sub, double T1, double T2, int n)
{
	int i;
	double diff, Tcur, Cv;
	
	diff = (T2-T1)/((double)n-1.0);
	for(i=0; i<n; i++) {
		Tcur = T1+diff*i;
		Cv = cvEval(sub, Tcur);
		printf("%g %g\n", Tcur, Cv);
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
	double q1[3], q2[3],T1,T2, *ret, Cv;
	//Index ordering is coupled with subInd()
	const sp subs[4] = {Ne, Ar, Kr, Xe};
	if(argc < 2)
		error("Cmon man... This waay to few arguments, ARGUE MORE!\n");
		
	if((subindex = subInd(argv[1]))<0)
		error("Unkown (for the program) substance\n");
	
	if((mode = operMode(argv[2]))<0)
		error("Mode is not implemented\n");
	
	//Dealing with the command line arguments
	switch(mode) {
		case 0: //Omega
			switch(argc) {
				case 6:
					//We only give frequency for one point
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					ret = freqEval(subs[subindex], q1);
					printVal(q1,ret);
					free(ret);
					break;
				case 9:
					//Two points, frequency need to be evaluated in 11 points
					//Equally spaced
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nEval(subs[subindex], q1, q2, 11, &freqEval);
					break;
				case 10:
					//Two points, frequency needs to evaluated in <npoints> points
					//Equally spaced between q1 and q2
					npoints = atof(argv[9]);
					if(npoints<2)
						error("Need to evaluate two or more points if a range is given\n");
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nEval(subs[subindex], q1, q2, npoints, &freqEval);
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
					ret = volDepEval(subs[subindex], q1);
					if(ret != NULL)
						printVal(q1,ret);
					else
						error("Not defined at [0 0 0]\n");
					break;
				case 9:
					//Two points, need to be spaced 11
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nEval(subs[subindex], q1, q2, 11, &volDepEval);
					break;
				case 10:
					//Two points, volDep needs to evaluated in <npoints> points
					//Equally spaced between q1 and q2
					npoints = atof(argv[9]);
					if(npoints<2)
						error("Need to evaluate two or more points if a range is given\n");
					q1[0] = atof(argv[3]); q1[1] = atof(argv[4]); q1[2] = atof(argv[5]);
					q2[0] = atof(argv[6]); q2[1] = atof(argv[7]); q2[2] = atof(argv[8]);
					nEval(subs[subindex], q1, q2, npoints, &volDepEval);
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		case 2: // cv
			switch(argc) {
				case 4:
					//One point calculate Cv
					T1 = atof(argv[3]);
					Cv = cvEval(subs[subindex],T1);
					printf("%g %g\n", T1, Cv);
					break;
				case 5:
					//Two points, need to be spaced 11
					T1 = atof(argv[3]);
					T2 = atof(argv[4]);
					nCvEval(subs[subindex], T1, T2, 11);
					break;
				case 6:
					//Two Temperatures needs to evaluated in <npoints> points
					//Equally spaced between T1 and T2
					npoints = atof(argv[5]);
					if(npoints<2)
						error("Need to evaluate two or more points if a range is given\n");
					T1 = atof(argv[3]);
					T2 = atof(argv[4]);
					nCvEval(subs[subindex], T1, T2, npoints); 
					break;
				default:
					error("Wrong number of arguments\n");
			}
			break;
		default:
			error("Congratz, you've reached what the programmer thought to be an unreachable part of the program\n");
	}
	
	return 0;	
}
