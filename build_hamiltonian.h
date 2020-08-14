#pragma once

#include <thread>

#include "parameters.h"
#include "ham_functions.h"
#include "matrix.h"
#include "basis.h"

using std::thread;

Matrix RotationAndTrap(vector<int>& channels, unsigned int QN);
Matrix isoInteraction(vector<int>& channels, unsigned int QN, double iso_strength);
Matrix anisoInteraction(vector<int>& channels, unsigned int QN, double aniso_strength);
Matrix electricField(vector<int>& channels, unsigned int QN, double el_field);
Matrix magneticField(vector<int>& channels, double mg_field);
Matrix spinRotation(vector<int>& channels, double spinrot_strength);
void addToHamiltonian(Matrix& hamiltonian, Matrix& matrix);
void buildHamiltonian(Matrix& hamiltonian, vector<int>& channels, unsigned int QN, double gjj, double gjk, double E, double B, double A, bool rotation_and_trap = true);
Matrix buildHamiltonianFast(vector<int>& channels, unsigned int QN, double gjj, double gjk, double E, double B, double A);
double getElectricFieldCoefficient(int j1, int j1prim, int j2, int j2prim, int J, int Jprim, int M, int Mprim, double el_field);
double getInteractionCoefficient(int n, int nprim, double interaction);
double getSpinRotationCoefficient(int j1, int j2, int J, int M, int Stot, int MS, int j1prim, int j2prim, int Jprim, int Mprim, int Stotprim, int MSprim, double spinrot_strength);