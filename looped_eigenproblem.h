#pragma once

#include "parameters.h"
#include "eigenproblem.h"
#include "utility.h"
#include "build_hamiltonian.h"

void loopedEigenproblem(Matrix& base_hamiltonian, vector<int>& channels, unsigned int QN, vector<string>& params, unsigned int how_many_eigenstates_to_save);
void findEnergiesAndSave(double value, Matrix& matrix_to_add, int iterator, vector<string>& params, Matrix& hamiltonian, Matrix& base_hamiltonian, unsigned int basis_size, unsigned int how_many_eigenstates_to_save, bool print = false);
void printHeigenvectors(pair <double*, double*>& eigenproblem, vector<int>& channels, unsigned int QN, vector<string>& params, int iterator = 0);
vector <double> findExpectedJ2(vector <pair <double, int> >& set_of_pairs, pair <double*, double*>& eigenproblem, vector<int>& channels, unsigned int QN);
vector <double> findExpectedMg(vector <pair <double, int> >& set_of_pairs, pair <double*, double*>& eigenproblem, vector<int>& channels, unsigned int QN, bool print = false);
vector < pair <double, int> > pairAndSortEigenenergies(pair <double*, double*>& eigenproblem, unsigned int basis_size);
int indexOfSpinBasisState(int n, int J, int M, int j1, int j2, int Stot, int MS, vector<int>& channels, unsigned int QN);
void printBasisState(int n, int J, int M, int j1, int j2, int Stot = 5, int MS = 5);
void printBasisState(int index, vector<int>& channels, unsigned int QN);
void printHamiltonianElement(int i, int j, vector<int>& channels, Matrix& hamiltonian, unsigned int QN);
void printChosenEigenvector(int state, pair <double*, double*>& eigenproblem, vector<int>& channels, unsigned int QN, bool sorted = true, std::ostream& ostr = cout);