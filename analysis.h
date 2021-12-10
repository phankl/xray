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

vector<double> moments(vector<double>);

vector<vector<double>> projectionMatrix(int);
vector<double> backwardSubstitution(vector<vector<double>>, vector<double>);

void save(vector<double>, vector<double>, string);

#endif
