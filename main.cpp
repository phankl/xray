#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

#include "analysis.h"
#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
 
  int nTube = 1;
  vector<double> tubeProperties = {0.6785e1, 100.0};
  vector<int> rasterProperties = {50, 1000};
  XYZ boxSize(20000.0, 20000.0, 20000.0);
  XYZ voxelSize(20000, 20000, 20000);
  XYZ mean(0.0, 0.0, 1000.0);
  XYZ std(0.001, 0.001, 0.001);

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

  int nBins = 200;
  int nMoments = 30;
  double angleSpacing2D = 2.0 * numbers::pi / nBins;
  double angleSpacing3D = numbers::pi / nBins;
  vector<double> angles2D(nBins, 0.0);
  vector<double> angles3D(nBins, 0.0);
  for (int i = 0; i < nBins; i++) {
    angles2D[i] = (i+0.5) * angleSpacing2D;
    angles3D[i] = (i+0.5) * angleSpacing3D;
  }

  // converging distribution

  int nStructure = 1;
  vector<double> intensityDistribution2D(nBins, 0.0);
  vector<double> odf3D(nBins, 0.0);
  vector<vector<double>> intensity2D(voxelSize.x, vector<double>(voxelSize.z, 0.0));
 
  for (int i = 0; i < nStructure; i++) {
    Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, i, true);
    vector<double> intensityTemp = structure.intensityDistribution2D(qRadius, nBins);
    vector<double> odfTemp = structure.odf(nBins);
    for (int j = 0; j < nBins; j++) {
      intensityDistribution2D[j] += intensityTemp[j] / nStructure;
      odf3D[j] += odfTemp[j] / nStructure;
    }
    for (int j = 0; j < voxelSize.x; j++)
      for (int k = 0; k < voxelSize.z; k++)
        intensity2D[j][k] += structure.intensity2D[j][k] / nStructure;
    cout << "Generated structure number " << i+1 << endl;
  }
  
  Structure structure(nTube, rasterProperties, tubeProperties, boxSize, voxelSize, mean, std, 1, true);
  //structure.intensity2D = intensity2D;
  //structure.saveIntensity("intensity2D.vtk");
  save(angles2D, intensityDistribution2D, "intensityDistribution2D.dat");
  save(angles3D, odf3D, "odf.dat");

  // moment calculation and inverse projection
  
  vector<double> intensityHalf2D(nBins/2, 0.0);
  for (int i = 0; i < nBins/2; i++)
    intensityHalf2D[i] = intensityDistribution2D[i];
  vector<double> intensityNaiveMoments3D = moments3D(intensityHalf2D, nMoments);

  vector<vector<double>> projection = projectionMatrix(nMoments);
  vector<double> intensityMoments2D = moments2D(intensityDistribution2D, nMoments);
  vector<double> intensityMoments3D = backwardSubstitution(projection, intensityMoments2D);
  
  vector<double> intensityDistributionTube3D = structure.intensityDistributionTube(qRadius, nBins);
  vector<double> intensityMomentsTube3D = moments3D(intensityDistributionTube3D, nMoments);
  vector<double> odfMoments3D = moments3D(odf3D, nMoments);

  save(angles3D, intensityDistributionTube3D, "intensityDistributionTube3D.dat");

  vector<double> odfMomentsReconstructed3D(nMoments, 0.0);

  for (int i = 0; i < nMoments; i += 2) {
    odfMomentsReconstructed3D[i] = intensityMoments3D[i] / intensityMomentsTube3D[i];
    double ratio = odfMomentsReconstructed3D[i] / odfMoments3D[i];
    cout << i << ": " << intensityMoments2D[i] << " " 
                      << intensityNaiveMoments3D[i] << " "
                      << intensityMoments3D[i] << " " 
                      << intensityMomentsTube3D[i] << " " 
                      << odfMomentsReconstructed3D[i] << " " 
                      << odfMoments3D[i] << " " 
                      << ratio << endl;
  }
  
  vector<double> intensityDistributionReconstructed2D(nBins, 0.0);
  for (int i = 0; i < nBins; i++) {
    double angle = (i+0.5) * angleSpacing2D;
    intensityDistributionReconstructed2D[i] += 0.5*intensityMoments2D[0] / numbers::pi;
    for (int j = 1; j < nMoments; j++)
      intensityDistributionReconstructed2D[i] += intensityMoments2D[j] * cos(j * angle) / numbers::pi;
  }

  save(angles2D, intensityDistributionReconstructed2D, "intensityDistributionReconstructed2D.dat");

  return 0;
}
