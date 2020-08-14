#pragma once

#include "mkl_lapacke.h"
#include <utility>
#include <iostream>

#include "constants.h"
#include "matrix.h"

using std::make_pair;
using std::pair;
using std::cout;

pair <double*, double*> MKLDiagonalization(Matrix& matrix, unsigned int matrix_size, bool compute_eigenvectors = false);
extern void printEigenvalues(const char* desc, MKL_INT n, double* w);
extern void printEigenvectors(const char* desc, MKL_INT n, double* v);