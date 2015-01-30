#include <stdio.h>
#include <gsl/gsl_matrix.h>
int main(){
 gsl_matrix * M;
 M = gsl_matrix_calloc (10, 12);
 gsl_matrix_set (M, 1, 3, 5.);
 gsl_matrix_fprintf (stdout, M, "%e");
 return 0;
}
