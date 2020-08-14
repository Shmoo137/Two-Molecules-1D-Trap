#pragma once

// in trap frequency units

const double PI = 3.14159265358979323846;
const double gs = 2.001159652181;
const double B1 = PI;
const double B2 = PI;
const double mu1 = 1;
const double mu2 = 1;

const unsigned int QN_nospin = 5; // number of quantum numbers, here: n, Jtot, Mtot, j1, j2
const unsigned int QN_withspin = 9; // number of quantum numbers, here: n, Jtot, Mtot, j1, j2, s1, ms1, s2, ms2
const int NMAX = 20;
const int JMAX = 4;
const int S1 = 1;
const int S2 = S1;
const double dS1 = (double)S1 * 0.5;
const double dS2 = (double)S2 * 0.5;

const bool ENERGY_DIFFERENCES = false; // instead of spectrum, the energy differences are calculated. Useful for quench analysis.