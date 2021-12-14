#include "analysis.h"

long double binomial(int n, int k) {
  if (k == 0 || n == k) return 1;
  else return binomial(n-1, k-1) * n / k; 
}

double sinc(double x) {
  if (x == 0.0) return 1.0;
  return sin(x) / x;
}

vector<double> moments2D(vector<double> intensity, int nMoments) {

  /*
  // calculate fourier coefficients of intensity distribution using fft

  // setup fft plan

  int n = intensity.size();
  vector<double> input(n, 0.0);
  vector<complex<double>> resultComplex(n, {0.0, 0.0});
  fftw_complex *fftwResult = reinterpret_cast<fftw_complex*>(&resultComplex[0]);

  fftw_plan plan = fftw_plan_dft_r2c_1d(n, &input[0], fftwResult, FFTW_ESTIMATE);
  input = intensity;
  
  // run fft and clean up

  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  // calculate moments
  
  vector<double> result(nMoments, 0.0);
  for (int i = 0; i < nMoments; i++) result[i] = 2.0 * numbers::pi * resultComplex[i].real() / n;
   for (int i = nMoments-1; i >= 0; i--) result[i] /= result[0];
  */
  
  vector<double> result(nMoments, 0.0);
  double angleSpacing = 2.0 * numbers::pi / intensity.size();

  for (int i = 0; i < nMoments; i++) {
    double sum = 0.0;
    for (int j = 0; j < intensity.size(); j++) {
      double angle = (j+0.5) * angleSpacing;
      sum += intensity[j] * cos(i * angle);
    }
    sum *= angleSpacing;
    result[i] = sum;
  }

  for (int i = nMoments-1; i >= 0; i--)
    result[i] /= result[0];

  return result;
}

vector<double> moments3D(vector<double> intensity, int nMoments) {

  vector<double> result(nMoments, 0.0);
  double angleSpacing = numbers::pi / intensity.size();

  for (int i = 0; i < nMoments; i++) {
    double sum = 0.0;
    for (int j = 0; j < intensity.size(); j++) {
      double angle = (j+0.5) * angleSpacing;
      sum += intensity[j] * sin(angle) * legendre(i, cos(angle));
    }
    sum *= angleSpacing;
    result[i] = sum;
  }

  for (int i = nMoments-1; i >= 0; i--)
    result[i] /= result[0];

  return result;
}

vector<vector<double>> projectionMatrix(int n) {
  
  vector<vector<double>> result(n, vector<double>(n, 0.0));
  result[0][0] = 1.0;
  
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j += 2) {
      long double prefactor = 0.5 * sqrt(numbers::pi) * (2*j + 1) / pow(4, j);
      long double sum = 0.0;
      for (int k = 0; k <= (j-i)/2; k++) {
        long double temp = pow(-4, k) * binomial(j, k) * binomial(2*j-2*k, j) * binomial(j-2*k, (j-i-2*k)/2);
        temp *= tgamma(0.5*(j-2*k+2)) / tgamma(0.5*(j-2*k+3));
        sum += temp;
      }
      result[i][j] = prefactor * sum;
    }

  return result;
}

vector<double> backwardSubstitution(vector<vector<double>> a, vector<double> b) {
  
  int n = a.size();
  vector<double> x(n, 0.0);
  
  x[n-1] = b[n-1] / a[n-1][n-1];

  for (int i = n-2; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i+1; j < n; j++)
      sum += a[i][j] * x[j];
    x[i] = (b[i] - sum) / a[i][i];
  }

  return x;
}

void save(vector<double> x, vector<double> y, string fileName) {
  ofstream file(fileName);

  for (int i = 0; i < x.size(); i++)
    file << x[i] << " " << y[i] << endl;
  
  file.close();
}
