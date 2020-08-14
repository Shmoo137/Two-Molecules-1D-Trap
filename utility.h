#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <iomanip> // setprecision
#include <algorithm> // sort
#include <fstream>
#include <math.h>

#include "constants.h"
#include "parameters.h"
#include "basis.h"

using std::vector;
using std::cout;
using std::setprecision;
using std::endl;
using std::string;
using std::to_string;
using std::pair;

void readArguments(int argc, char* argv[]);
void saveResults(string what_quantity, double* w, unsigned int basis_size, vector<string>& params, int iterator, double value, unsigned int how_many);
void saveResults(string what_quantity, vector<double>& results, vector<string>& params, int iterator, double value, unsigned int how_many);
vector < pair <double, int> > pairAndSort(vector<double>& first_set);
vector < pair <int, int> > pairAndSort(vector <pair <double, int> >& set_of_pairs_second);
vector < pair <double, int> > pairAndSort(vector<double>& first_set, vector<int>& second_set);
vector < pair <double, double> > pairAndSort(vector<double>& first_set, vector<double>& second_set);
string createFilename(vector<string>& params, int iterator, string what_quantity);
string createQuenchFilename(vector<string>& params, vector<string>& quenched_params, string what_quantity, bool time = false);
vector <string> argToStrings(char* argv[]);
vector <string> quenchedToString();
string doubleToString(double number);
void coutSquareMatrix(vector<double>& matrix, int dim, int precision);
void saveSquareMatrix(vector<double>& matrix, int dim, int precision, string filename);
void saveVector(vector<double>& vector, int precision, string filename);