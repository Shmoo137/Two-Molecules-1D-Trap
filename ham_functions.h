#pragma once

#include "constants.h"
#include <cmath>
#include <stdlib.h>
#include <vector>
using std::vector; //since this is a header file, we don’t need to use the whole namespace

long double factorial(unsigned int liczba, unsigned int m);
double Hermite0(int n);
long double ClebschGordanCoefficient(int j1, int m1, int j2, int m2, int Jtot, int Mtot);
long double SpinClebschGordanCoefficient(double j1, double m1, double j2, double m2, int Jtot, int Mtot);
double diagonalFunction(int n, int j1, int j2);
double electricFunction(double mu, double el_field, int j, int jprim, int m);
double magneticFunction(double mg_field, double sum_of_mss);
double deltaFunction(double g, int n1, int n2);
double spinRotationDiagonalFunction(double spinrot_strength, int m1, int m2, double ms1, double ms2);

bool differByNOnly(vector<int>& channels, int QN, int i, int j); // for isoInteraction
bool differByNAndjsOnly(vector<int>& channels, int QN, int i, int j); // for anisoInteraction
bool jsExchangeByOne(vector<int>& channels, int QN, int i, int j); // for anisoInteraction
bool sameNAndSpin(vector<int>& channels, int QN, int i, int j); // for electricField
bool onejDifferByOne(vector<int>& channels, int QN, int i, int j); // for electricField
bool differByJMjSAndMsOnly(vector<int>& channels, int i, int j); // for spinRotation
bool msAndmssOffDiagonalSpinRotation(int m1, int m1prim, int m2, int m2prim, double ms1, double ms1prim, double ms2, double ms2prim); // for spinRotation
bool nonZeroCG(int m1, int m2, int M, int m1prim, int m2prim, int Mprim); // for electricField and spinRotation
bool nonZeroSpinCG(double m1, double m2, double M, double m1prim, double m2prim, double Mprim); // for spinRotation
bool momentumConserved(vector<int>& channels, int QN, int i, int j); // for anisoInteraction