#include "analysis.h"

void save(vector<double> x, vector<double> y, string fileName) {
  ofstream file(fileName);

  for (int i = 0; i < x.size(); i++)
    file << x[i] << " " << y[i] << endl;
  
  file.close();
}
