#include "structure.h"

Structure::Structure(int nTubeNew, vector<int> rasterProperties, vector<double> tubeProperties, XYZ boxSizeNew, XYZ nVoxelNew, XYZ meanNew, XYZ stdNew):
  nTube(nTubeNew),
  rTube(tubeProperties[0]),
  lTube(tubeProperties[1]),
  phiPoints(rasterProperties[0]),
  zPoints(rasterProperties[1]),
  boxSize(boxSizeNew),
  nVoxel(nVoxelNew),
  mean(meanNew),
  std(stdNew)
{
  generate();
  voxelise();
  fft();
}

void Structure::savePoints(string fileName) {
  ofstream file(fileName);

  file << points.size() << endl;
  file << "Atoms. Timestep: 0" << endl;

  for (int i = 0; i < points.size(); i++) {
    XYZ point = points[i];
    file << 1 << " " << point.x << " " << point.y << " " << point.z << endl;
  }
  file.close();
}

void Structure::saveDensity(string fileName) {
  ofstream file(fileName);

  // VTK file header

  file << "# vtk DataFile Version 4.2" << endl;
  file << "Charge density voxel data" << endl;
  file << "ASCII" << endl;

  // define data structure

  file << "DATASET STRUCTURED_POINTS" << endl;
  file << "DIMENSIONS " << nVoxel.x << " " << nVoxel.y << " " << nVoxel.z << endl;
  file << "ORIGIN 0.0 0.0 0.0" << endl;
  file << "SPACING " << voxelSize.x << " " << voxelSize.y << " " << voxelSize.z << endl;

  // data block

  file << "POINT_DATA " << int(nVoxel.x * nVoxel.y * nVoxel.z) << endl;
  file << "SCALARS DENSITY FLOAT 1" << endl;
  file << "LOOKUP_TABLE DENSITY" << endl;

  for (int k = 0; k < nVoxel.z; k++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int i = 0; i < nVoxel.z; i++)
        file << density[i][j][k] << endl;

  file.close();
}

void Structure::saveIntensity(string fileName) {
  ofstream file(fileName);

  // VTK file header

  file << "# vtk DataFile Version 4.2" << endl;
  file << "X-ray intensity voxel data" << endl;
  file << "ASCII" << endl;

  // define data structure

  file << "DATASET STRUCTURED_POINTS" << endl;
  file << "DIMENSIONS " << nVoxel.x << " " << nVoxel.y << " " << nVoxel.z << endl;
  file << "ORIGIN 0.0 0.0 0.0" << endl;
  file << "SPACING " << voxelSize.x << " " << voxelSize.y << " " << voxelSize.z << endl;

  // data block

  file << "POINT_DATA " << int(nVoxel.x * nVoxel.y * nVoxel.z) << endl;
  file << "SCALARS INTENSITY FLOAT 1" << endl;
  file << "LOOKUP_TABLE INTENSITY" << endl;

  for (int k = 0; k < nVoxel.z; k++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int i = 0; i < nVoxel.z; i++)
        file << intensity[i][j][k] << endl;

  file.close();
}

double Structure::distance(Tube tube1, Tube tube2, double lambda1, double lambda2) {
  return (tube1.s - tube2.s + lambda1*tube1.t - lambda2*tube2.t).length();
}

bool Structure::checkOverlap(Tube tube1) {
  
  double r1 = tube1.r;
  double l1 = tube1.l;

  XYZ t1 = tube1.t;
  XYZ s1 = tube1.s;

  for (int i = 0; i < tubes.size(); i++) {   
    
    Tube tube2 = tubes[i];
    
    double r2 = tube2.r;
    double l2 = tube2.l;

    XYZ t2 = tube2.t;
    XYZ s2 = tube2.s;

    // check for minimum distance of lines

    XYZ n = cross(t1, t2);
    n.normalise();

    XYZ delta = s2 - s1;

    double d = fabs(delta * n);

    if (d > r1+r2) continue;

    // check if minimum distance is within line segments

    double c = t1 * t2;
    double frac = 1.0 / (1.0 - c*c);
    
    double lambda1 = (t1 - c*t2) * delta * frac;
    double lambda2 = (c*t1 - t2) * delta * frac;
 
    if (lambda1 >= 0.0 && lambda1 <= l1 && lambda2 >= 0.0 && lambda2 <= l2) {
      if (d < r1+r2) return true;
      else continue;
    }

    // check cases where minimum distance is outside the line segments
    // not fully exploiting convexity on subdomains yet, kinda brute force

    // domain edges
    
    lambda1 = 0.0;
    lambda2 = -(delta * t2);
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = c*l1 - delta*t2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = delta * t1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    lambda1 = c*l2 + delta*t1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    // domain vertices

    lambda1 = 0.0;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = 0.0;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
    
    lambda1 = l1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
  }

  return false;
}

void Structure::generate() {

  mt19937_64 generator(2021);
  generator.discard(10000);

  uniform_real_distribution<double> sDistribution(0.0, 1.0);
  normal_distribution<double> tDistributionX(mean.x, std.x);
  normal_distribution<double> tDistributionY(mean.y, std.y);
  normal_distribution<double> tDistributionZ(mean.z, std.z);

  int attempts = 0;
  int tubeIndex = 1;

  while (tubes.size() < nTube) {
    
    // generate new tube

    double sX = sDistribution(generator) * boxSize.x;
    double sY = sDistribution(generator) * boxSize.y;
    double sZ = sDistribution(generator) * boxSize.z;
  
    double tX = tDistributionX(generator);
    double tY = tDistributionY(generator);
    double tZ = tDistributionZ(generator);
  
    XYZ s(sX, sY, sZ);
    XYZ t(tX, tY, tZ);

    t.normalise();

    Tube tube(rTube, lTube, s, t);

    attempts++;

    // add new tube if no overlap is detected

    if (!checkOverlap(tube)) {
      cout << "Tube " << tubeIndex << " generated after " << attempts << " attempts." << endl;
      attempts = 0;
      tubes.push_back(tube);
      tubeIndex++;
    }
  }
}

void Structure::voxelise() {
    
  double phiSpacing = 6.28318530718 / phiPoints;
  
  // generate points for charge density

  for (int i = 0; i < nTube; i++) {
    
    Tube tube = tubes[i];
    XYZ s = tube.s;
    XYZ t = tube.t;
    double r = tube.r;
    double l = tube.l;

    double zSpacing = l / (zPoints - 1);

    // generate points along surface of tube

    XYZ n;
    if (t.x != 0.0) n = XYZ(-(t.y + t.z)/t.x, 1.0, 1.0);
    else if (t.y != 0.0) n = XYZ(1.0, -(t.x + t.z)/t.y, 1.0);
    else n = XYZ(1.0, 1.0, -(t.x + t.y)/t.z);
    n.normalise();

    for (int j = 0; j < phiPoints; j++) {

      XYZ nRotated = rotate(n, t, j*phiSpacing);
      nRotated.normalise();
      
      for (int k = 0; k < zPoints; k++) {
        XYZ point = s + r*nRotated + k*zSpacing*t;
        
        // apply periodic boundary conditions
        
        if (point.x < 0.0) point.x += boxSize.x;
        if (point.y < 0.0) point.y += boxSize.y;
        if (point.z < 0.0) point.z += boxSize.z;
        if (point.x > boxSize.x) point.x -= boxSize.x;
        if (point.y > boxSize.y) point.y -= boxSize.y;
        if (point.z > boxSize.z) point.z -= boxSize.z;

        points.push_back(point);
      } 
    }
  }

  // translate points to voxelised charge density

  double voxelSizeX = boxSize.x / nVoxel.x;
  double voxelSizeY = boxSize.y / nVoxel.y;
  double voxelSizeZ = boxSize.z / nVoxel.z;

  voxelSize = XYZ(voxelSizeX, voxelSizeY, voxelSizeZ);

  density = vector<vector<vector<double>>>(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));

  for (int i = 0; i < points.size(); i++) {
    XYZ point = points[i];

    int indexX = floor(point.x / voxelSizeX);
    int indexY = floor(point.y / voxelSizeY);
    int indexZ = floor(point.z / voxelSizeZ);
  
    density[indexX][indexY][indexZ] += 1.0;
  }
}

void Structure::fft() {
  
  // setup fft procedure
  
  vector<double> densityUnwrapped(int(nVoxel.x*nVoxel.y*nVoxel.z), 0.0);
  vector<complex<double>> amplitudeUnwrapped(int(nVoxel.x*nVoxel.y)*(int(nVoxel.z)/2 + 1), {0.0, 0.0});
  fftw_complex *fftwAmplitudeUnwrapped = reinterpret_cast<fftw_complex*>(&amplitudeUnwrapped[0]);

  fftw_plan plan = fftw_plan_dft_r2c_3d(int(nVoxel.x), int(nVoxel.y), int(nVoxel.z), &densityUnwrapped[0], fftwAmplitudeUnwrapped, FFTW_ESTIMATE);

  // unwrap density data to 1D vector

  for (int i = 0; i < nVoxel.x; i++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int k = 0; k < nVoxel.z; k++)
        densityUnwrapped[i*int(nVoxel.y)*int(nVoxel.z) + j*int(nVoxel.z) + k] = density[i][j][k];

  // run fft and clean up

  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // calculate intensity and wrap up into 3d vector

  intensity = vector<vector<vector<double>>>(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));
  
  int nVoxelZTemp = int(nVoxel.z/2) + 1;

  for (int i = 0; i < nVoxel.x; i++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int k = 0; k < nVoxel.z; k++) {
        int kTemp = k;
        if (k > int(nVoxel.z)/2) kTemp = int(nVoxel.z) - k - 1;
        int index = i*int(nVoxel.y)*nVoxelZTemp + j*nVoxelZTemp + kTemp;
        intensity[i][j][k] = norm(amplitudeUnwrapped[index]);
      } 

  // wrap fft around and center at origin

  
  vector<vector<vector<double>>> intensityCentred(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));
  
  for (int i = 0; i < nVoxel.x; i++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int k = 0; k < nVoxel.y; k++) {
        int iTemp, jTemp, kTemp;
        if (i > nVoxel.x/2) iTemp = i - int(nVoxel.x)/2;
        else iTemp = i + int(nVoxel.x)/2 - 1;
        if (j > nVoxel.y/2) jTemp = j - int(nVoxel.y)/2;
        else jTemp = j + int(nVoxel.y)/2 - 1;
        if (k > nVoxel.z/2) kTemp = k - int(nVoxel.z)/2;
        else kTemp = k + int(nVoxel.z)/2 - 1;
      
        intensityCentred[i][j][k] = intensity[iTemp][jTemp][kTemp];
      }

  intensity = intensityCentred;
}
