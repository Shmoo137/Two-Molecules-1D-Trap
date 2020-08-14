#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
using std::vector; //since this is a header file, we don’t need to use the whole namespace
using std::pair;

#include "constants.h"
#include "parameters.h"
#include "ham_functions.h"
#include "basis.h"

double r2MeanValue(int n, int nprim);
double j_r2_jprim(int j, int jprim, pair <double*, double*>& eigenproblem, vector<int>& channels, int QN);
double j_Jtot2_jprim(int j, int jprim, pair <double*, double*>& eigenproblem, vector<int>& channels, int QN);
double j_j2_jprim(int j, int jprim, pair <double*, double*>& eigenproblem, vector<int>& channels, int QN);
double j_Mg_jprim(int j, int jprim, pair <double*, double*>& eigenproblem, vector<int>& channels, int QN);
double j_Observable_jprim(int j, int jprim, pair <double*, double*>& eigenproblem, vector<int>& channels, int QN);
bool indexOccursInVector(unsigned int index, vector<unsigned int>& vector_of_indices);