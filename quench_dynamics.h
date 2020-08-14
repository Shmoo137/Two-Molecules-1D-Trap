#pragma once

#include "parameters.h"
#include "quench_functions.h"
#include "matrix.h"
#include "looped_eigenproblem.h"
#include "utility.h"

vector <double> chosenInitialStateBasisCoefficients(Matrix& hamiltonian_initial, vector<int>& channels, unsigned int QN, int which_state = 0);
vector <double> calculatePhiPsiJOverlaps(pair <double*, double*>& eigenproblem_final, vector<double>& chosen_initial_state_coefficients, unsigned int basis_size, vector<string>& string_params, vector<string>& string_quenched_params, bool save_overlaps = false);
vector <double> calculateMatrixOfPsiOPsiOverlaps(pair <double*, double*>& eigenproblem_final, vector<int>& channels, unsigned int QN);
vector <unsigned int> findConvergedPsiJIndices(pair <double*, double*>& eigenproblem_final, vector<int>& channels, unsigned int QN, bool print = false);
void quenchDynamics(Matrix& hamiltonian_final, vector<int>& channels, unsigned int QN, vector<double>& chosen_initial_state_coefficients, vector<string>& string_params, vector<string>& string_quenched_params, bool save_couplings = true);