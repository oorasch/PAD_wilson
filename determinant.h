#ifndef FILE_DETERMINANT
#define FILE_DETERMINANT
#include <vector>
#include <cmath>
#include <algorithm>
#include<gsl/gsl_linalg.h>

/*
 *
 * name: det
 * @param 
 * @return
 *
 */
double det(const int& n, std::vector<double>& mat);

/*
 *
 * name: det
 * @param 
 * @return
 *
 */
double ln_det(const int& n, std::vector<double>& mat);


#endif
