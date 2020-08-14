#pragma once // non-standard but widely supported preprocessor directive designed to cause the current source file to be included only once in a single compilation
// built with this tutorial: https://medium.com/@furkanicus/how-to-create-a-matrix-class-using-c-3641f37809c7

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <string>

using std::vector; //since this is a header file, we don’t need to use the whole namespace
using std::string; //since this is a header file, we don’t need to use the whole namespace

class Matrix {
	private: // it ensures that these members can be used only inside the class methods.They will be inaccessible from main.cpp
		unsigned m_rowSize;
		unsigned m_colSize;

	protected:
		vector<vector<double> > m_matrix;

	public: // constructors - they initiate the code inside the main.cpp
		Matrix(unsigned, unsigned, double); // take the row size, column size and an initial value for each cell

		// Matrix Operations
		// we're overloading operators +, - and *
		Matrix operator+(Matrix&); // Matrix indicates that output will be of Matrix type, Matrix& - that the operation is done on the matrix itself, not its copy
		Matrix operator-(Matrix&);
		Matrix operator*(Matrix&);
		void operator+=(Matrix&);
		Matrix transpose();

		// Scalar Operations
		Matrix operator+(double);
		Matrix operator-(double);
		Matrix operator*(double);
		Matrix operator/(double);

		// Aesthetic Methods
		double& operator()(const unsigned&, const unsigned&); // This line enables to write it in this form: Matrix(3,4)
		double& operator()(const unsigned&); // This line enables to get a diagonal argument (3,3) in the form: Matrix(3)
		void print(string text = "Matrix: ") const; // print the whole matrix in Matlab form. Takes an optional argument "text header" and outputs nothing
		unsigned getRows() const;
		unsigned getCols() const;
};