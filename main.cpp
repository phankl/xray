#include <iostream>
#include <fstream>
#include <vector>

#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
  
  int nTube = 1000;
  vector<double> tubeProperties = {1.0, 100.0};
  vector<int> rasterProperties = {20, 400};
  XYZ boxSize(200.0, 200.0, 200.0);
  XYZ voxelSize(500, 500, 500);
  XYZ mean(0.0, 0.0, 1.0);
  XYZ std(1.0, 1.0, 3.0);

  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std);
  structure.savePoints("points.xyz");
  structure.saveDensity("density.vtk");
  structure.saveIntensity("intensity.vtk");

  return 0;
}
