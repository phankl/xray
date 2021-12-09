#ifndef ANALYSIS_HEADER
#define ANALYSIS_HEADER

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

vector<double> moments(vector<double>, int);

vector<vector<double>> projectionMatrix(int);
vector<vector<double>> forwardSubstitution(vector<vector<double>>);

void save(vector<double>, vector<double>, string);

#endif
