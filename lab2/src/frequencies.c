#include <math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

void
frequencies(double A, double B, double m, double *q,
				double *omega, double *eps)
{
	int i,j;
	double val;
	gsl_matrix* D;
	gsl_eigen_symmv_workspace* w;
	gsl_vector* ovec;
	gsl_matrix* emat;
	
	//We only consider 3x3 matrices
	D = gsl_matrix_alloc(3,3);
	
	//I really don't care about efficiency here
	for(i=0; i<3; i++) {
		//printf("%d %d %d\n", i, (i+1)%3, (i+2)%3);
		gsl_matrix_set(D, i, i, 0.5/m*(A+B)*(8-4*cos(q[i]*M_PI)*cos(q[(i+1)%3]*M_PI)
						- 4*cos(q[i]*M_PI)*cos(q[(i+2)%3]*M_PI))
						+ B/m*(4-4*cos(q[(i+1)%3]*M_PI)*cos(q[(i+2)%3]*M_PI)));
		
		gsl_matrix_set(D, i, (i+1)%3, 0.5/m*(A-B)*4*sin(q[i]*M_PI)*sin(q[(i+1)%3]*M_PI));
		gsl_matrix_set(D, i, (i+2)%3, 0.5/m*(A-B)*4*sin(q[i]*M_PI)*sin(q[(i+2)%3]*M_PI));
	}
	
	//Solving
	ovec = gsl_vector_alloc(3);
	emat = gsl_matrix_alloc(3,3);
	w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv(D, ovec, emat, w);

	//sort
	gsl_eigen_symmv_sort(ovec,emat,GSL_EIGEN_SORT_VAL_ASC);
	
	//Translate values
	for(i=0; i<3; i++) {
		omega[i] = fmax(0.0,sqrt(gsl_vector_get(ovec,i)));
		if(eps != NULL)
			for(j=0; j<3; j++)
				eps[i+j*3] = gsl_matrix_get(emat,i,j);
	}

	//free
	gsl_matrix_free(D);gsl_matrix_free(emat);
	gsl_vector_free(ovec);
	gsl_eigen_symmv_free(w);
}
