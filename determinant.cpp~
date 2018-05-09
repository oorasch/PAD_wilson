#include "determinant.h"

using namespace std;

double det(const int& n, vector<double>& mat)
{
	int s;
	gsl_matrix_view m = gsl_matrix_view_array (mat.data(), n, n);
	gsl_permutation * p = gsl_permutation_alloc (n);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  
  gsl_permutation_free(p);

 	return gsl_linalg_LU_det(&m.matrix, s);

}


double ln_det(const int& n, vector<double>& mat)
{
	int s;
	gsl_matrix_view m = gsl_matrix_view_array (mat.data(), n, n);
	gsl_permutation * p = gsl_permutation_alloc (n);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  
  gsl_permutation_free(p);

 	return gsl_linalg_LU_lndet(&m.matrix);
}
