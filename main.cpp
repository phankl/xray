#include <iostream>
#include <fstream>
#include <vector>

#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
  
  int nTube = 3000;
  vector<double> tubeProperties = {1.0, 50.0};
  vector<int> rasterProperties = {20, 400};
  XYZ boxSize(200.0, 200.0, 200.0);
  XYZ voxelSize(5000, 5000, 5000);
  XYZ mean(0.0, 0.0, 5.0);
  XYZ std(1.0, 1.0, 1.0);

  // SAXS true

  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, true);
  //structure.savePoints("points.xyz");
  structure.saveDensity("density2D.vtk");
  structure.saveIntensity("intensity2D.vtk");

  // SAXS false
  /*  
  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, false);
  // structure.savePoints("points.xyz");
  // structure.saveDensity("density.vtk");
  structure.saveIntensity("intensity.vtk");
  structure.ewaldSphere(XYZ(300.0, 0.0, 0.0), "ewaldSphere.vtk");
  */

  return 0;
}
