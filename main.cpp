#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

#include "analysis.h"
#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
  
  int nTube = 10000;
  vector<double> tubeProperties = {0.6785, 100.0};
  vector<int> rasterProperties = {20, 400};
  XYZ boxSize(1000.0, 1000.0, 1000.0);
  XYZ voxelSize(1000, 1000, 1000);
  XYZ mean(0.0, 0.0, 5.0);
  XYZ std(1.0, 1.0, 1.0);

  // SAXS true

  //Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, true);
  //structure.savePoints("points.xyz");
  //structure.saveDensity("density2D.vtk");
  //structure.saveIntensity("intensity2D.vtk");

  // SAXS false
  /*  
  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, false);
  // structure.savePoints("points.xyz");
  // structure.saveDensity("density.vtk");
  structure.saveIntensity("intensity.vtk");
  structure.ewaldSphere(XYZ(300.0, 0.0, 0.0), "ewaldSphere.vtk");
  */

  // output analysis

  double qRadius = 0.5;

  int nBins = 100;
  double angleSpacing2D = 2.0 * numbers::pi / nBins;
  vector<double> angles2D(nBins, 0.0);
  for (int i = 0; i < nBins; i++) angles2D[i] = (i+0.5) * angleSpacing2D;
  
  // converging distribution

  int nStructure = 100;
  vector<double> intensityDistribution2D(nBins, 0.0);
  vector<vector<double>> intensity2D(voxelSize.x, vector<double>(voxelSize.z, 0.0));
 
  for (int i = 0; i < nStructure; i++) {
    Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, i, true);
    vector<double> intensityTemp = structure.intensityDistribution2D(qRadius, nBins);
    for (int j = 0; j < nBins; j++) intensityDistribution2D[j] += intensityTemp[j] / nStructure;
    for (int j = 0; j < voxelSize.x; j++)
      for (int k = 0; k < voxelSize.z; k++)
        intensity2D[j][k] += structure.intensity2D[j][k] / nStructure;
    cout << "Generated structure number " << i+1 << endl;
  }
  
  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, 1, true);
  structure.intensity2D = intensity2D;
  structure.saveIntensity("intensity2D.vtk");
  save(angles2D, intensityDistribution2D, "intensityDistribution2D.dat");

  // moment calculation and inverse projection

  vector<double> moments2D = moments(intensityDistribution2D);
  vector<vector<double>> projection = projectionMatrix(nBins/2 + 1);
  vector<double> moments3D = backwardSubstitution(projection, moments2D);

  return 0;
}
