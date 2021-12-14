#ifndef ANALYSIS_HEADER
#define ANALYSIS_HEADER

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <numbers>

#include <complex>
#include <fftw3.h>

using namespace std;

long double binomial(int, int);
double sinc(double);

vector<double> moments2D(vector<double>, int);
vector<double> moments3D(vector<double>, int);

vector<vector<double>> projectionMatrix(int);
vector<double> backwardSubstitution(vector<vector<double>>, vector<double>);

void save(vector<double>, vector<double>, string);

#endif
